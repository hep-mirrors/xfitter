
/*
   @file ReactionBaseDISNC.cc
   @date 2017-04-08
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-04-08
*/

#include "ReactionBaseDISNC.h"
#include <iostream>
#include <cstdio>
#include "QCDNUM/QCDNUM.h"
#include "QCDNUM_Manager.h"
#include <IntegrateDIS.h>
#include "xfitter_pars.h"
#include "hf_errlog.h"
#include "BaseEvolution.h"
#include "EvolutionQCDNUM.h"
#include <spline.h>

template <typename T>
void print(T d)
{
  std::cout << d << "\n";
}

// Helpers for QCDNUM:

//! F2,FL full
const double CNEP2F[] = {0., 0., 1., 0., 1., 0., 0., 0., 1., 0., 1., 0., 0.}; //u  (top off ?)
const double CNEM2F[] = {0., 1., 0., 1., 0., 1., 0., 1., 0., 1., 0., 1., 0.}; //d

//! xF3 full
const double CNEP3F[] = {0., 0., -1., 0., -1., 0., 0., 0., 1., 0., 1., 0., 0.};  //u
const double CNEM3F[] = {0., -1., 0., -1., 0., -1., 0., 1., 0., 1., 0., 1., 0.}; //d

//! c
const double CNEP2Fc[] = {0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.}; //c

//! b
const double CNEM2Fb[] = {0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.}; //b

// define QCDNUM function:
extern "C"
{
  void zmstfun_(const int &id, const double &key, double &x, double &q2, double &sf, const int &np, const int &flag);
}

// the class factories
extern "C" ReactionBaseDISNC *create()
{
  return new ReactionBaseDISNC();
}

// Initialize at the start of the computation
void ReactionBaseDISNC::atStart()
{
}

// Main function to compute results at an iteration
void ReactionBaseDISNC::compute(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal)
{
  unsigned termID = td->id;

  valarray<double> val;
  map<string, valarray<double>> err;

  switch (GetDataType(termID))
  {
  case dataType::signonred:
  {
    sred(td, val, err);
    // transform reduced -> non-reduced cross sections
    auto &x = *GetBinValues(td, "x");
    auto &q2 = *GetBinValues(td, "Q2");
    auto &y = *GetBinValues(td, "y");
    const double pi = 3.1415926535897932384626433832795029;
    valarray<double> yplus = 1.0 + (1.0 - y) * (1.0 - y);
    valarray<double> factor = 2 * pi * _alphaem * _alphaem * yplus / (q2 * q2 * x) * _convfac;
    val *= factor;
    break;
  }
  case dataType::sigred:
    sred(td, val, err);
    break;
  case dataType::f2:
    F2(td, val, err);
    if (_ht[termID]) 
      ApplyHigherTwist(td, 2, val, err);
    break;
  case dataType::fl:
    FL(td, val, err);
    if (_ht[termID]) 
      ApplyHigherTwist(td, 1, val, err);
    break;
  }

  if (_integrated.find(termID) == _integrated.end())
  {
    // usual cross section at (q2,x) points
    valExternal = val;
    errExternal = err;
  }
  else
  {
    // integrated cross sections
    valExternal = _integrated[termID]->compute(val);
    // no idea how error could be treated: for now do nothing
    errExternal = err;
  }
}

void ReactionBaseDISNC::atIteration()
{
  // Make sure to call the parent class initialization:
  super::atIteration();

  _convfac = *XFITTER_PARS::getParamD("convFac");
  _alphaem = *XFITTER_PARS::getParamD("alphaem");
  _Mz = *XFITTER_PARS::getParamD("Mz");
  _Mw = *XFITTER_PARS::getParamD("Mw");
  _sin2thetaW = *XFITTER_PARS::getParamD("sin2thW");

  _ve = -0.5 + 2. * _sin2thetaW; // !
  _ae = -0.5;                    // !
  _au = 0.5;
  _ad = -0.5;
  _vu = _au - (4. / 3.) * _sin2thetaW;
  _vd = _ad + (2. / 3.) * _sin2thetaW;

  //  print (_Mz);

  // Re-set internal maps (faster access):
  for (auto ds : _dsIDs)
  {
    (_f2u[ds])[0] = -100.;
    (_flu[ds])[0] = -100.;
    (_xf3u[ds])[0] = -100.;
  }
}

//
void ReactionBaseDISNC::initTerm(TermData *td)
{
  unsigned termID = td->id;

  {
  const string& name = getReactionName();
  if (name == "BaseDISNC" or name == "RT_DISNC") {
    //Reaction requires QCDNUM, make sure it is used
    xfitter::BaseEvolution* pdf = td->getPDF();
    if ( pdf->getClassName() != string("QCDNUM") ){
      cerr<<"[ERROR] "<<getReactionName()<<" can only work with QCDNUM evolution; got evolution \""<<pdf->_name<<"\" of class \""<<pdf->getClassName()<<"\" for termID="<<termID<<endl;
      hf_errlog(19052311, "F: Chosen DISNC reaction can only work with QCDNUM evolution, see stderr for details");
    }

    xfitter::requireZMSTF();
  }
  }

  _dsIDs.push_back(termID);

  _tdDS[termID] = td;
  _polarisation[termID] = (td->hasParam("epolarity")) ? *td->getParamD("epolarity") : 0;
  _charge[termID] = (td->hasParam("echarge")) ? *td->getParamD("echarge") : 0;

  // Inclusive reduced cross section by default.
  _dataType[termID] = dataType::sigred;
  _dataFlav[termID] = dataFlav::incl;
  string msg = "I: Calculating DIS NC reduced cross section";

  if (td->hasParam("type"))
  {
    string type = td->getParamS("type");
    if (type == "signonred")
    {
      _dataType[termID] = dataType::signonred;
      msg = "I: Calculating DIS NC double-differential (non-reduced) cross section";
    }
    else if (type == "sigred")
    {
      _dataType[termID] = dataType::sigred;
      msg = "I: Calculating DIS NC reduced cross section";
    }
    else if (type == "F2")
    {
      _dataType[termID] = dataType::f2;
      msg = "I: Calculating DIS NC F2";
    }
    else if (type == "FL")
    {
      _dataType[termID] = dataType::fl;
      msg = "I: Calculating DIS NC FL";
    }
    else
    {
      char buffer[256];
      sprintf(buffer, "F: dataset with id = %d has unknown type = %s", termID, type.c_str());
      string str = buffer;
      hf_errlog(17101901, str);
    }
  }

  // flav: incl, c, b
  if (td->hasParam("flav"))
  {
    string flav = td->getParamS("flav");
    if (flav == "incl")
    {
      _dataFlav[termID] = dataFlav::incl;
      msg += " inclusive";
    }
    else if (flav == "c")
    {
      _dataFlav[termID] = dataFlav::c;
      msg += " charm";
    }
    else if (flav == "b")
    {
      _dataFlav[termID] = dataFlav::b;
      msg += " beauty";
    }
    else
    {
      char buffer[256];
      sprintf(buffer, "F: dataset with id = %d has unknown flav = %s", termID, flav.c_str());
      string str = buffer;
      hf_errlog(17101902, str);
    }
  }

  // check if centre-of-mass energy is provided
  double s = -1.0;
  if (td->hasParam("energy"))
    s = pow(*td->getParamD("energy"), 2.0);

  // bins
  // if Q2min, Q2max, ymin and ymax (and optionally xmin, xmax) are provided, calculate integrated cross sections
  auto *q2minp = td->getBinColumnOrNull("Q2min");
  auto *q2maxp = td->getBinColumnOrNull("Q2max");
  // also try small first letter for Q2 (for backward compatibility)
  if (!q2minp)
    q2minp = td->getBinColumnOrNull("q2min");
  if (!q2maxp)
    q2maxp = td->getBinColumnOrNull("q2max");
  auto *yminp = td->getBinColumnOrNull("ymin");
  auto *ymaxp = td->getBinColumnOrNull("ymax");
  // optional xmin, xmax for integrated cross sections
  auto *xminp = td->getBinColumnOrNull("xmin");
  auto *xmaxp = td->getBinColumnOrNull("xmax");

  if (q2minp && q2maxp)
  {
    // integrated cross section
    if (s < 0)
      hf_errlog(18060100, "F: centre-of-mass energy is required for integrated DIS dataset " + std::to_string(termID));
    if (_dataType[termID] != dataType::signonred)
      hf_errlog(18060200, "F: integrated DIS can be calculated only for non-reduced cross sections, dataset " + std::to_string(termID));
    IntegrateDIS *iDIS = new IntegrateDIS();
    _npoints[termID] = iDIS->init(s, q2minp, q2maxp, yminp, ymaxp, xminp, xmaxp);
    _integrated[termID] = iDIS;
    msg += " (integrated)";
  }
  else
  {
    // cross section at (Q2,x) points
    auto *q2p = td->getBinColumnOrNull("Q2"), *xp = td->getBinColumnOrNull("x"), *yp = td->getBinColumnOrNull("y");

    // if Q2 and x bins and centre-of-mass energy provided, calculate y = Q2 / (s * x)

    /*
          DO NOT ALLOW FOR MISSING Y COLUMN

    if(yp == nullptr && q2p != nullptr && xp != nullptr)
    {
      if ( s > 0.0 )
      {
        valarray<double> y = (*q2p) / (s * (*xp));
        std::pair<string,valarray<double>* > dsBin = std::make_pair("y", &y);
        AddBinning(termID, &dsBin);
        yp = GetBinValues(termID, "y");
      }
    }
    */

    if (q2p == nullptr || xp == nullptr || yp == nullptr)
    {
      string msg = "F: Q2, x or y bins are missing for NC DIS reaction for dataset " + std::to_string(termID);
      hf_errlog(17040801, msg);
    }
    _npoints[termID] = (*q2p).size();
  }
  if (td->hasParam("ht")) {
    _ht[termID] = td->getParamI("ht");
  }

  hf_errlog(17041001, msg);

  // Allocate internal arrays:
  _f2u[termID].resize(_npoints[termID]);
  _f2d[termID].resize(_npoints[termID]);
  _flu[termID].resize(_npoints[termID]);
  _fld[termID].resize(_npoints[termID]);
  _xf3u[termID].resize(_npoints[termID]);
  _xf3d[termID].resize(_npoints[termID]);

  // Get PDF id
  _ipdfSet[termID] = static_cast<xfitter::EvolutionQCDNUM*> (td->getPDF())->getPdfType();

  // higher twist spline knots
  _ht_x = {    0.,    0.1,    0.3,   0.5,   0.7,   0.9, 1.};
  _ht_2 = { 0.023, -0.032, -0.005, 0.025, 0.051, 0.003, 0.};
  _ht_t = {-0.319, -0.134, -0.052, 0.071, 0.030, 0.003, 0.};
  _ht_alpha_2 = 0.;
  _ht_alpha_t = 0.05;
}

void ReactionBaseDISNC::reinitTerm(TermData *td)
{
  unsigned termID = td->id;
  _polarisation[termID] = (td->hasParam("epolarity")) ? *td->getParamD("epolarity") : 0;
}

const valarray<double> *ReactionBaseDISNC::GetBinValues(TermData *td, const string &binName)
{
  unsigned termID = td->id;

  if (_integrated.find(termID) == _integrated.end())
    return td->getBinColumnOrNull(binName);
  else
  {
    if (binName == "Q2")
      return _integrated[termID]->getBinValuesQ2();
    else if (binName == "x")
      return _integrated[termID]->getBinValuesX();
    else if (binName == "y")
      return _integrated[termID]->getBinValuesY();
    else
      return td->getBinColumnOrNull(binName);
  }
};

void ReactionBaseDISNC::F2gamma BASE_PARS
{
  valarray<double> f2u, f2d;
  GetF2ud(td, f2u, f2d);
  val = 4. / 9. * f2u + 1. / 9. * f2d;
}

void ReactionBaseDISNC::F2gammaZ BASE_PARS
{
  valarray<double> f2u, f2d;
  GetF2ud(td, f2u, f2d);
  val = 2. * (2. / 3. * _vu * f2u - 1. / 3. * _vd * f2d);
}

void ReactionBaseDISNC::F2Z BASE_PARS
{
  valarray<double> f2u, f2d;
  GetF2ud(td, f2u, f2d);
  val = (_vu * _vu + _au * _au) * f2u + (_vd * _vd + _ad * _ad) * f2d;
}

void ReactionBaseDISNC::F2 BASE_PARS
{
  unsigned termID = td->id;

  valarray<double> f2g(_npoints[termID]);
  F2gamma(td, f2g, err);

  valarray<double> f2gZ(_npoints[termID]);
  F2gammaZ(td, f2gZ, err);

  valarray<double> f2Z(_npoints[termID]);
  F2Z(td, f2Z, err);

  valarray<double> k(_npoints[termID]);
  kappa(td, k);
  // combine together:

  double pol = GetPolarisation(termID);
  double charge = GetCharge(termID);

  val = f2g - (_ve + charge * pol * _ae) * k * f2gZ + (_ae * _ae + _ve * _ve + 2 * charge * pol * _ae * _ve) * k * k * f2Z;
}

void ReactionBaseDISNC::FLgamma BASE_PARS
{
  valarray<double> flu, fld;
  GetFLud(td, flu, fld);
  val = 4. / 9. * flu + 1. / 9. * fld;
}

void ReactionBaseDISNC::FLgammaZ BASE_PARS
{
  valarray<double> flu, fld;
  GetFLud(td, flu, fld);
  val = 2. * (2. / 3. * _vu * flu - 1. / 3. * _vd * fld);
}

void ReactionBaseDISNC::FLZ BASE_PARS
{
  valarray<double> flu, fld;
  GetFLud(td, flu, fld);
  val = (_vu * _vu + _au * _au) * flu + (_vd * _vd + _ad * _ad) * fld;
}

void ReactionBaseDISNC::FL BASE_PARS
{
  unsigned termID = td->id;

  valarray<double> flg(_npoints[termID]);
  FLgamma(td, flg, err);

  valarray<double> flgZ(_npoints[termID]);
  FLgammaZ(td, flgZ, err);

  valarray<double> flZ(_npoints[termID]);
  FLZ(td, flZ, err);

  valarray<double> k(_npoints[termID]);
  kappa(td, k);
  // combine together:

  double pol = GetPolarisation(termID);
  double charge = GetCharge(termID);

  val = flg - (_ve + charge * pol * _ae) * k * flgZ + (_ae * _ae + _ve * _ve + 2 * charge * pol * _ae * _ve) * k * k * flZ;
}

void ReactionBaseDISNC::xF3gammaZ BASE_PARS
{
  valarray<double> xf3u, xf3d;
  GetxF3ud(td, xf3u, xf3d);
  val = 2. * (2. / 3. * _au * xf3u - 1. / 3. * _ad * xf3d);
}

void ReactionBaseDISNC::xF3Z BASE_PARS
{
  valarray<double> xf3u, xf3d;
  GetxF3ud(td, xf3u, xf3d);
  val = 2. * (_vu * _au * xf3u + _vd * _ad * xf3d);
}

void ReactionBaseDISNC::xF3 BASE_PARS
{
  unsigned termID = td->id;
  valarray<double> xf3gZ(_npoints[termID]);
  xF3gammaZ(td, xf3gZ, err);

  valarray<double> xf3Z(_npoints[termID]);
  xF3Z(td, xf3Z, err);

  valarray<double> k(_npoints[termID]);
  kappa(td, k);

  double pol = GetPolarisation(termID);
  double charge = GetCharge(termID);

  val = (_ae * charge + pol * _ve) * k * xf3gZ + (-2 * _ae * _ve * charge - pol * (_ve * _ve + _ae * _ae)) * k * k * xf3Z;
}

void ReactionBaseDISNC::sred BASE_PARS
{
  unsigned termID = td->id;

  auto *yp = GetBinValues(td, "y");
  auto y = *yp;

  valarray<double> f2(_npoints[termID]);
  F2(td, f2, err);
  if (_ht[termID]) 
    ApplyHigherTwist(td, 2, f2, err);

  valarray<double> fl(_npoints[termID]);
  FL(td, fl, err);
  if (_ht[termID]) 
    ApplyHigherTwist(td, 1, fl, err);

  valarray<double> xf3(_npoints[termID]);
  xF3(td, xf3, err);

  //  double charge = GetCharge(termID);   xF3 is alredy charge-dependent.

  valarray<double> yplus = 1.0 + (1.0 - y) * (1.0 - y);
  valarray<double> yminus = 1.0 - (1.0 - y) * (1.0 - y);

  val = f2 - y * y / yplus * fl + (yminus / yplus) * xf3;
}

void ReactionBaseDISNC::GetF2ud(TermData *td, valarray<double> &f2u, valarray<double> &f2d)
{
  // Check if already computed:
  unsigned termID = td->id;

  if ((_f2u[termID])[0] < -99.)
  { // compute
    // Get x,Q2 arrays:
    auto *q2p = GetBinValues(td, "Q2"), *xp = GetBinValues(td, "x");
    auto q2 = *q2p, x = *xp;

    // Call QCDNUM
    QCDNUM::zswitch(_ipdfSet[termID]);       // This sets proper PDF set for the computations below

    const int id = 2;
    const int flag = 0;
    int Npnt = GetNpoint(termID);
    switch (GetDataFlav(termID))
    {
    case dataFlav::incl:
      zmstfun_(id, CNEP2F[0], x[0], q2[0], (_f2u[termID])[0], Npnt, flag);
      zmstfun_(id, CNEM2F[0], x[0], q2[0], (_f2d[termID])[0], Npnt, flag);
      break;
    case dataFlav::c:
      zmstfun_(id, CNEP2Fc[0], x[0], q2[0], (_f2u[termID])[0], Npnt, flag);
      break;
    case dataFlav::b:
      zmstfun_(id, CNEM2Fb[0], x[0], q2[0], (_f2d[termID])[0], Npnt, flag);
      break;
    }
  }
  f2u = _f2u[termID];
  f2d = _f2d[termID];
}

void ReactionBaseDISNC::GetFLud(TermData *td, valarray<double> &flu, valarray<double> &fld)
{
  // Check if already computed:
  unsigned termID = td->id;
  if ((_flu[termID])[0] < -99.)
  { // compute
    // Get x,Q2 arrays:
    auto *q2p = GetBinValues(td, "Q2"), *xp = GetBinValues(td, "x");
    auto q2 = *q2p, x = *xp;

    // Call QCDNUM
    QCDNUM::zswitch(_ipdfSet[termID]);       // This sets proper PDF set for the computations below
    const int id = 1;
    const int flag = 0;
    int Npnt = GetNpoint(termID);
    switch (GetDataFlav(termID))
    {
    case dataFlav::incl:
      zmstfun_(id, CNEP2F[0], x[0], q2[0], (_flu[termID])[0], Npnt, flag);
      zmstfun_(id, CNEM2F[0], x[0], q2[0], (_fld[termID])[0], Npnt, flag);
      break;
    case dataFlav::c:
      zmstfun_(id, CNEP2Fc[0], x[0], q2[0], (_flu[termID])[0], Npnt, flag);
      break;
    case dataFlav::b:
      zmstfun_(id, CNEM2Fb[0], x[0], q2[0], (_fld[termID])[0], Npnt, flag);
      break;
    }
  }
  flu = _flu[termID];
  fld = _fld[termID];
}

void ReactionBaseDISNC::GetxF3ud(TermData *td, valarray<double> &xf3u, valarray<double> &xf3d)
{
  // Check if already computed:
  unsigned termID = td->id;
  if ((_xf3u[termID])[0] < -99.)
  { // compute
    // Get x,Q2 arrays:
    auto *q2p = GetBinValues(td, "Q2"), *xp = GetBinValues(td, "x");
    auto q2 = *q2p, x = *xp;

    // Call QCDNUM
    QCDNUM::zswitch(_ipdfSet[termID]);       // This sets proper PDF set for the computations below
    const int id = 3;
    const int flag = 0;
    int Npnt = GetNpoint(termID);
    // OZ 19.10.2017 TODO: F3 is 0 in VFNS for heavy quarks?
    //if ( GetDataType(termID) == dataType::sigred ) {
    if (GetDataFlav(termID) == dataFlav::incl)
    {
      zmstfun_(id, CNEP3F[0], x[0], q2[0], (_xf3u[termID])[0], Npnt, flag);
      zmstfun_(id, CNEM3F[0], x[0], q2[0], (_xf3d[termID])[0], Npnt, flag);
    }
    else
    {
      for (int p = 0; p < Npnt; p++)
        _xf3u[termID][p] = _xf3d[termID][p] = 0;
    }
    //else {
    //  NOT_IMPLEMENTED("xF3 c,b");
    //}
  }
  xf3u = _xf3u[termID];
  xf3d = _xf3d[termID];
}

void ReactionBaseDISNC::kappa(TermData *td, valarray<double> &k)
{
  auto *q2p = GetBinValues(td, "Q2");
  double cos2thetaW = 1 - _sin2thetaW;

  k = 1. / (4 * _sin2thetaW * cos2thetaW) * (*q2p) / ((*q2p) + _Mz * _Mz);
}

void ReactionBaseDISNC::ApplyHigherTwist(TermData *td, const int f_type, valarray<double>& val, map<string, valarray<double>>& err)
{
  if (f_type == 1) {
    // F_t = F_2 - F_L
    valarray<double> f2;
    F2(td, f2, err);
    valarray<double> ft = f2 - val;
    tk::spline spline;
    spline.set_points(_ht_x, _ht_t);
    auto &x = *GetBinValues(td, "x");
    auto &q2 = *GetBinValues(td, "Q2");
    for (size_t ip = 0; ip < ft.size(); ip++) {
      ft[ip] += pow(x[ip], _ht_alpha_t) * spline(x[ip]) / q2[ip];
    }
    // F_L = F_2 - F_T
    val = f2 - ft;
  }
  else if (f_type == 2) {
    tk::spline spline;
    spline.set_points(_ht_x, _ht_2);
    auto &x = *GetBinValues(td, "x");
    auto &q2 = *GetBinValues(td, "Q2");
    for (size_t ip = 0; ip < val.size(); ip++) {
      val[ip] += pow(x[ip], _ht_alpha_2) * spline(x[ip]) / q2[ip];
    }
  }
}
