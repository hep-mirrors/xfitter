
/*
   @file ReactionBaseDISCC.cc
   @date 2017-10-05
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-10-05
*/

#include "ReactionBaseDISCC.h"
#include <iostream>
#include "QCDNUM/QCDNUM.h"
#include "QCDNUM_Manager.h"
#include "xfitter_pars.h"
#include "hf_errlog.h"
#include "BaseEvolution.h"
#include "EvolutionQCDNUM.h"

// Helpers for QCDNUM (CC):

//! full
const double CCEP2F[] = {0., 0., 1., 0., 1., 0., 0., 1., 0., 1., 0., 0., 0.};
const double CCEM2F[] = {0., 0., 0., 1., 0., 1., 0., 0., 1., 0., 1., 0., 0.};

const double CCEP3F[] = {0., 0., -1., 0., -1., 0., 0., 1., 0., 1., 0., 0., 0.};
const double CCEM3F[] = {0., 0., 0., -1., 0., -1., 0., 0., 1., 0., 1., 0., 0.};

//! c
// work in progress: according to 1001.2312 section 5,
// in ZM only the sum of contributions s + c makes sense
// three different options are below for checks, uncommented one is for s + c
//
// only c
//const double  CCEP2Fc[] = {0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.} ;
//const double  CCEM2Fc[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.} ;
// only s
//const double  CCEP2Fc[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.} ;
//const double  CCEM2Fc[] = {0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.} ;
// only s,c
const double CCEP2Fc[] = {0., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.};
const double CCEM2Fc[] = {0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 0.};

// only c
//const double  CCEP3Fc[] = {0.,0.,-1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
//const double  CCEM3Fc[] = {0.,0. ,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.};
// only s
//const double  CCEP3Fc[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.};
//const double  CCEM3Fc[] = {0.,0.,0.,-1.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
// only s,c
const double CCEP3Fc[] = {0., 0., -1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.};
const double CCEM3Fc[] = {0., 0., 0., -1., 0., 0., 0., 0., 0., 0., 1., 0., 0.};

// define QCDNUM function:
extern "C"
{
//key, x, q2, sf are arrays, better use pointer rather than reference
void zmstfun_(const int &id, const double &key, double &x, double &q2, double &sf, const int &np, const int &flag);
}

// the class factories
extern "C" ReactionBaseDISCC *create()
{
  return new ReactionBaseDISCC();
}

// Initialize at the start of the computation
void ReactionBaseDISCC::atStart()
{
}

valarray<double> GetF(TermData *td, const int id)
{
  /*F   id incl   c       isNegative
  FL  1  CCEP2F CCEP2Fc F
  F2  2  CCEP2F CCEP2Fc F
  xF3 3  CCEP3F CCEP3Fc F
  FL  1  CCEM2F CCEM2Fc T
  F2  2  CCEM2F CCEM2Fc T
  xF3 3  CCEM3F CCEM3Fc T
*/
  BaseDISCC::ReactionData *rd = (BaseDISCC::ReactionData *)td->reactionData;
  bool isNegative = rd->_charge < 0;
  auto &q2 = *BaseDISCC::GetBinValues(td, "Q2"),
      &x = *BaseDISCC::GetBinValues(td, "x");
  const double *C;
  switch (rd->_dataFlav)
  {
    case BaseDISCC::dataFlav::incl:
      if (isNegative)
      {
        if (id == 3)
          C = CCEM3F;
        else
          C = CCEM2F;
      }
      else
      {
        if (id == 3)
          C = CCEP3F;
        else
          C = CCEP2F;
      }
      break;
    case BaseDISCC::dataFlav::c:
      if (isNegative)
      {
        if (id == 3)
          C = CCEM3Fc;
        else
          C = CCEM2Fc;
      }
      else
      {
        if (id == 3)
          C = CCEP3Fc;
        else
          C = CCEP2Fc;
      }
      break;
    case BaseDISCC::dataFlav::b:
      hf_errlog(2019060601,"F:Flavor b is not implemented in BaseDISCC");
    default:
      std::abort(); //unreachable
  }
  // Call QCDNUM
  QCDNUM::zswitch(rd->ipdfSet); // Set PDF in QCDNUM:
  const int flag = 0;
  const int Npnt = x.size();
  valarray<double> ret(Npnt);
  zmstfun_(id, C[0], const_cast<double &>(x[0]), const_cast<double &>(q2[0]), ret[0], Npnt, flag);
  return ret;
}

valarray<double> ReactionBaseDISCC::FL(TermData *td) { return GetF(td, 1); }
valarray<double> ReactionBaseDISCC::F2(TermData *td) { return GetF(td, 2); }
valarray<double> ReactionBaseDISCC::xF3(TermData *td) { return GetF(td, 3); }

// Main function to compute results at an iteration
void ReactionBaseDISCC::compute(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal)
{
  BaseDISCC::ReactionData *rd = (BaseDISCC::ReactionData *)td->reactionData;
  const double MW = *rd->Mw;

  // Basic formulae for CC cross section:
  const valarray<double> &y = *BaseDISCC::GetBinValues(td, "y");

  valarray<double> yplus = 1.0 + (1.0 - y) * (1.0 - y);
  valarray<double> yminus = 1.0 - (1.0 - y) * (1.0 - y);

  valarray<double> f2 = F2(td);
  valarray<double> fl = FL(td);
  valarray<double> xf3 = xF3(td);
  double polarity = rd->_polarisation;
  double charge = rd->_charge;

  valarray<double> val;
  xf3 = 0.;
  //f2 = fl = 0.;
  if (charge > 0)
    val = 0.5 * (1 + polarity) * (yplus * f2 - yminus * xf3 - y * y * fl);
  else
    val = 0.5 * (1 - polarity) * (yplus * f2 + yminus * xf3 - y * y * fl);

  if (!rd->_isReduced)
  {
    // extra factor for non-reduced cross section
    auto &x = *BaseDISCC::GetBinValues(td, "x"),
        &q2 = *BaseDISCC::GetBinValues(td, "Q2");
    const double pi = 3.1415926535897932384626433832795029;
    valarray<double> factor = (MW * MW * MW * MW / pow((q2 + MW * MW), 2)) * _Gf * _Gf / (2 * pi * x) * _convfac;
    val *= factor;
  }

  IntegrateDIS *iDIS = rd->_integrated;
  if (iDIS)
  {
    // integrated cross sections
    valExternal = iDIS->compute(val);
  }
  else
  {
    // usual cross section at (q2,x) points
    valExternal = val;
  }
}
void ReactionBaseDISCC::initTerm(TermData *td)
{
  unsigned termID = td->id;
  _dsIDs.push_back(termID);
  _tdDS[termID] = td;

  {
  const string& name = getReactionName();
  if (name == "BaseDISCC") {
    //Reaction requires QCDNUM, make sure it is used
    xfitter::BaseEvolution* pdf = td->getPDF();
    if ( pdf->getClassName() != string("QCDNUM") ){
      cerr<<"[ERROR] "<<getReactionName()<<" can only work with QCDNUM evolution; got evolution \""<<pdf->_name<<"\" of class \""<<pdf->getClassName()<<"\" for termID="<<termID<<endl;
      hf_errlog(19052312, "F: Chosen DISCC reaction can only work with QCDNUM evolution, see stderr for details");
    }

    xfitter::requireZMSTF();
  }
  }

  // This we do not want to fit:
  _Gf = *XFITTER_PARS::getParamD("gf");
  _convfac = *XFITTER_PARS::getParamD("convFac");

  BaseDISCC::ReactionData *rd = new BaseDISCC::ReactionData();
  td->reactionData = (void *)rd;
  auto &_isReduced = rd->_isReduced;
  auto &_dataFlav = rd->_dataFlav;
  auto& _npoints = rd->_npoints;
  if (td->hasParam("epolarity"))
    rd->_polarisation = *td->getParamD("epolarity"); //cannot be fitted
  if (td->hasParam("echarge"))
    rd->_charge = *td->getParamD("echarge"); //cannot be fitted
  if (td->hasParam("reduced"))
    _isReduced = td->getParamI("reduced");
  rd->Mw = td->getParamD("Mw");

  // Store QCDNUM evolution
  rd->ipdfSet = static_cast<xfitter::EvolutionQCDNUM*> (td->getPDF())->getPdfType();

  // type: sigred, signonred (no F2, FL implemented so far, thus type is defined by bool _isReduced)
  // HERA data files provide 'signonred' CC cross sections
  // Inclusive "non-reduced" cross section by default.
  string msg = "I: Calculating DIS CC reduced cross section";
  if (td->hasParam("type"))
  {
    string type = td->getParamS("type");
    if (type == "sigred")
    {
      _isReduced = true;
      msg = "I: Calculating DIS CC reduced cross section";
    }
    else if (type == "signonred")
    {
      _isReduced = false;
      msg = "I: Calculating DIS CC non-reduced cross section";
    }
    else
    {
      cerr << "[ERROR] Unknown type=\"" << type << "\" given to reaction \"" << getReactionName() << "\"; termID=" << td->id << endl;
      hf_errlog(17101903, "F: Unknown \"type\" given to reaction, see stderr");
    }
  }
  // flav: incl, c, b
  if (td->hasParam("flav"))
  {
    string flavor = td->getParamS("flav");
    if (flavor == "incl")
    {
      _dataFlav = BaseDISCC::dataFlav::incl;
      msg += " inclusive";
    }
    else if (flavor == "c")
    {
      _dataFlav = BaseDISCC::dataFlav::c;
      msg += " charm";
    }
    else if (flavor == "b")
    { //no beauty
      //NOT IMPLEMENTED
      hf_errlog(18042501, "F: predictions for beauty in CC are not available (term id = " + to_string(td->id) + ")");
    }
    else
    {
      cerr << "[ERROR] Unknown flavor=\"" << flavor << "\" given to reaction \"" << getReactionName() << "\"; termID=" << td->id << endl;
      hf_errlog(18042502, "F: Unknown \"flavor\" given to reaction, see stderr");
    }
  }
  // check if centre-of-mass energy is provided
  double s = -1.0;
  if (td->hasParam("energy"))
  {
    double energy = *td->getParamD("energy");
    s = energy * energy;
  }
  // bins
  // if Q2min, Q2max, ymin and ymax (and optionally xmin, xmax) are provided, integrated cross sections are calculated
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
      hf_errlog(18060100, "F: centre-of-mass energy is required for integrated DIS term " + std::to_string(td->id));
    if (_isReduced)
      hf_errlog(18060200, "F: integrated DIS can be calculated only for non-reduced cross sections, term " + std::to_string(td->id));
    IntegrateDIS *iDIS = new IntegrateDIS();
    _npoints = iDIS->init(s, q2minp, q2maxp, yminp, ymaxp, xminp, xmaxp);
    rd->_integrated = iDIS;
    msg += " (integrated)";
  }
  else
    _npoints = td->getNbins();
  hf_errlog(17041001, msg);
}

const valarray<double> *ReactionBaseDISCC::GetBinValues(TermData *td, const string &binName)
{
  //unsigned termID = td->id;
  BaseDISCC::ReactionData *rd = (BaseDISCC::ReactionData *)td->reactionData;

  if (!rd->_integrated)
    return td->getBinColumnOrNull(binName);
  else
  {
    if (binName == "Q2")
      return rd->_integrated->getBinValuesQ2();
    else if (binName == "x")
      return rd->_integrated->getBinValuesX();
    else if (binName == "y")
      return rd->_integrated->getBinValuesY();
    else
      return td->getBinColumnOrNull(binName);
  }
}
