 
/*
   @file ReactionBaseDISNC.cc
   @date 2017-04-08
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-04-08
*/

#include "ReactionBaseDISNC.h"
#include <iostream>
#include <cstdio>
#include <IntegrateDIS.h>

template <typename T>
void print(T d) {
  std::cout << d << "\n";
}

// Helpers for QCDNUM:

//! F2,FL full
const double CNEP2F[] = {0.,0.,1.,0.,1.,0.,0.,0.,1.,0.,1.,0.,0.}; //u  (top off ?)
const double CNEM2F[] = {0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0.}; //d

//! xF3 full
const double CNEP3F[] = {0., 0.,-1., 0.,-1., 0.,0.,0.,1.,0.,1.,0.,0.}; //u
const double CNEM3F[] = {0.,-1., 0.,-1., 0.,-1.,0.,1.,0.,1.,0.,1.,0.}; //d

//! c
const double CNEP2Fc[] = {0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.}; //c

//! b
const double CNEM2Fb[] = {0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.}; //b


// define QCDNUM function:
extern "C" {
  void zmstfun_(const int& id, const double& key, double& x, double& q2, double& sf, const int& np, const int &flag);
}


// the class factories
extern "C" ReactionBaseDISNC* create() {
  return new ReactionBaseDISNC();
}


// Initialize at the start of the computation
int ReactionBaseDISNC::initAtStart(const string &s)
{
  _convfac = GetParam("convFac");
  return 0;
}

// Main function to compute results at an iteration
int ReactionBaseDISNC::compute(int dataSetID, valarray<double> &valExternal, map<string, valarray<double> > &errExternal)
{
  valarray<double> val;
  map<string, valarray<double> > err;

  switch ( GetDataType(dataSetID) )
    {
    case dataType::signonred :
      {
        sred(dataSetID, val, err) ;
        // transform reduced -> non-reduced cross sections
        auto *xp  = GetBinValues(dataSetID,"x");
        auto x = *xp;
        auto *Q2p  = GetBinValues(dataSetID,"Q2");
        auto q2 = *Q2p;
        auto *yp  = GetBinValues(dataSetID,"y");
        auto y = *yp;
        const double pi = 3.1415926535897932384626433832795029;
        valarray<double> yplus  = 1.0+(1.0-y)*(1.0-y);
        valarray<double> factor = 2 * pi * _alphaem * _alphaem * yplus / (q2 * q2 * x) * _convfac;
        val *= factor;
        break ;
      }
    case dataType::sigred :
      sred(dataSetID, val, err) ;
      break ;
    case dataType::f2 :
      F2(dataSetID, val, err) ;
      break ;
    case dataType::fl :
      FL(dataSetID, val, err) ;
      break ;
    }

  if(_integrated.find(dataSetID) == _integrated.end())
  {
    // usual cross section at (q2,x) points
    valExternal = val;
    errExternal = err;
  }
  else
  {
    // integrated cross sections
    valExternal = _integrated[dataSetID]->compute(val);
    // no idea how error could be treated: for now do nothing
    errExternal = err;
  }

  return 0;
}

void ReactionBaseDISNC::initAtIteration() {
  _convfac = GetParam("convFac");
  _alphaem = GetParam("alphaem");
  _Mz = GetParam("Mz");
  _Mw = GetParam("Mw");
  _sin2thetaW = GetParam("sin2thW");
  
  _ve =  -0.5 + 2.*_sin2thetaW; // !
  _ae =  -0.5;                  // !
  _au =   0.5;
  _ad =  -0.5;
  _vu =  _au - (4./3.)*_sin2thetaW;
  _vd =  _ad + (2./3.)*_sin2thetaW;
  
  //  print (_Mz);

  // Re-set internal maps (faster access):
  for ( auto ds : _dsIDs)  {
    (_f2u[ds])[0] = -100.;
    (_flu[ds])[0] = -100.;
    (_xf3u[ds])[0] = -100.;
  }
}

// 
void  ReactionBaseDISNC::setDatasetParameters( int dataSetID, map<string,string> pars, map<string,double> parsDataset) 
{
  _polarisation[dataSetID] =  (parsDataset.find("epolarity") != parsDataset.end()) ? parsDataset["epolarity"] : 0;
  _charge[dataSetID]       =  (parsDataset.find("echarge")       != parsDataset.end()) ? parsDataset["echarge"] : 0;

  // Inclusive reduced cross section by default.
  _dataType[dataSetID] = dataType::sigred;
  _dataFlav[dataSetID] = dataFlav::incl;
  string msg = "I: Calculating DIS NC reduced cross section";
  if ( parsDataset.find("F2") != parsDataset.end() ) {
    _dataType[dataSetID] = dataType::f2;
    msg = "I: Calculating DIS NC F2";
  }
  if ( parsDataset.find("FL") != parsDataset.end() ) {
    _dataType[dataSetID] = dataType::fl;
    msg = "I: Calculating DIS NC FL";
  }
  if ( parsDataset.find("reduced") != parsDataset.end() ) {
    _dataType[dataSetID] = dataType::sigred;
    msg = "I: Calculating DIS NC reduced cross section";
  }

  // check if settings are provided in the new format key=value
  // type: signonred, sigred, F2, FL
  map<string,string>::iterator it = pars.find("type");
  if ( it != pars.end() ) {
    if(it->second == "signonred")
    {
      _dataType[dataSetID] = dataType::signonred;
      msg = "I: Calculating DIS NC double-differential (non-reduced) cross section";
    }
    else if(it->second == "sigred")
    {
      _dataType[dataSetID] = dataType::sigred;
      msg = "I: Calculating DIS NC reduced cross section";
    }
    else if(it->second == "F2")
    {
      _dataType[dataSetID] = dataType::f2;
      msg = "I: Calculating DIS NC F2";
    }
    else if(it->second == "FL")
    {
      _dataType[dataSetID] = dataType::fl;
      msg = "I: Calculating DIS NC FL";
    }
    else
    {
      char buffer[256];
      sprintf(buffer, "F: dataset with id = %d has unknown type = %s", dataSetID, it->second.c_str());
      string str = buffer;
      hf_errlog_(17101901, str.c_str(), str.length());
    }
  }

  // flav: incl, c, b
  it = pars.find("flav");
  if ( it != pars.end() ) {
    if(it->second == "incl")
    {
      _dataFlav[dataSetID] = dataFlav::incl;
      msg += " inclusive";
    }
    else if(it->second == "c")
    {
      _dataFlav[dataSetID] = dataFlav::c;
      msg += " charm";
    }
    else if(it->second == "b")
    {
      _dataFlav[dataSetID] = dataFlav::b;
      msg += " beauty";
    }
    else
    {
      char buffer[256];
      sprintf(buffer, "F: dataset with id = %d has unknown flav = %s", dataSetID, it->second.c_str());
      string str = buffer;
      hf_errlog_(17101902, str.c_str(), str.length());
    }
  }

  // e charge: double
  it = pars.find("echarge");
  if ( it != pars.end() )
    _charge[dataSetID] = atof(it->second.c_str());

  // e polarity: double
  it = pars.find("epolarity");
  if ( it != pars.end() )
    _polarisation[dataSetID] = atof(it->second.c_str());

  // check if centre-of-mass energy is provided
  double s = -1.0;
  map<string,string>::iterator itEnergy = pars.find("energy");
  if ( itEnergy != pars.end() )
    s = pow(stof(itEnergy->second), 2.0);

  // bins
  // if Q2min, Q2max, ymin and ymax (and optionally xmin, xmax) are provided, calculate integrated cross sections
  auto *q2minp  = GetBinValues(dataSetID,"Q2min");
  auto *q2maxp  = GetBinValues(dataSetID,"Q2max");
  // also try small first letter for backward compatibility
  if(!q2minp)
    q2minp  = GetBinValues(dataSetID,"q2min");
  if(!q2maxp)
    q2maxp  = GetBinValues(dataSetID,"q2max");
  auto *yminp  = GetBinValues(dataSetID,"ymin");
  auto *ymaxp  = GetBinValues(dataSetID,"ymax");
  // optional xmin, xmax for integrated cross sections
  auto *xminp  = GetBinValues(dataSetID,"xmin");
  auto *xmaxp  = GetBinValues(dataSetID,"xmax");

  if(q2minp && q2maxp && yminp && ymaxp)
  {
    // integrated cross section
    if(s < 0)
      hf_errlog(18060100, "F: centre-of-mass energy is required for integrated DIS dataset " + std::to_string(dataSetID));
    if(_dataType[dataSetID] != dataType::signonred)
      hf_errlog(18060200, "F: integrated DIS can be calculated only for non-reduced cross sections, dataset " + std::to_string(dataSetID));
    IntegrateDIS* iDIS = new IntegrateDIS();
    _npoints[dataSetID] = iDIS->init(s, q2minp, q2maxp, yminp, ymaxp, xminp, xmaxp);
    _integrated[dataSetID] = iDIS;
    msg += " (integrated)";
  }
  else
  {
    // cross section at (Q2,x) points
    auto *q2p  = GetBinValues(dataSetID,"Q2"), *xp  = GetBinValues(dataSetID,"x"), *yp  = GetBinValues(dataSetID,"y");

    // if Q2 and x bins and centre-of-mass energy provided, calculate y = Q2 / (s * x)
    if(yp == nullptr && q2p != nullptr && xp != nullptr)
    {
      if ( s > 0.0 )
      {
        valarray<double> y = (*q2p) / (s * (*xp));
        std::pair<string,valarray<double>* > dsBin = std::make_pair("y", &y);
        AddBinning(dataSetID, &dsBin);
        yp = GetBinValues(dataSetID, "y");
      }
    }

    if (q2p == nullptr || xp == nullptr || yp == nullptr ) {
      string msg = "F: Q2, x or Y bins are missing for NC DIS reaction for dataset " + std::to_string(dataSetID);
      hf_errlog_(17040801,msg.c_str(), msg.size());
    }
    _npoints[dataSetID] = (*q2p).size();
  }

  hf_errlog_(17041001, msg.c_str(), msg.size());

  // Allocate internal arrays:
  _f2u[dataSetID].resize(_npoints[dataSetID]);
  _f2d[dataSetID].resize(_npoints[dataSetID]);
  _flu[dataSetID].resize(_npoints[dataSetID]);
  _fld[dataSetID].resize(_npoints[dataSetID]);
  _xf3u[dataSetID].resize(_npoints[dataSetID]);
  _xf3d[dataSetID].resize(_npoints[dataSetID]);
}

valarray<double> *ReactionBaseDISNC::GetBinValues(int idDS, const string& binName)
{
  if(_integrated.find(idDS) == _integrated.end())
    return ReactionTheory::GetBinValues(idDS, binName);
  else
  {
    if(binName == "Q2")
      return _integrated[idDS]->getBinValuesQ2();
    else if(binName == "x")
      return _integrated[idDS]->getBinValuesX();
    else if(binName == "y")
      return _integrated[idDS]->getBinValuesY();
    else
      return ReactionTheory::GetBinValues(idDS, binName);
  }
};

void ReactionBaseDISNC::F2gamma BASE_PARS
{
  valarray<double> f2u, f2d;
  GetF2ud(dataSetID, f2u,  f2d);
  val = 4./9. * f2u + 1./9. * f2d;
}

void ReactionBaseDISNC::F2gammaZ BASE_PARS
{
  valarray<double> f2u, f2d;
  GetF2ud(dataSetID, f2u,  f2d);
  val = 2.*( 2./3.*_vu * f2u - 1./3.*_vd * f2d );
}

void ReactionBaseDISNC::F2Z BASE_PARS
{
  valarray<double> f2u, f2d;
  GetF2ud(dataSetID, f2u,  f2d);
  val = (_vu*_vu + _au*_au) * f2u + (_vd*_vd + _ad*_ad) * f2d ;
}

void ReactionBaseDISNC::F2 BASE_PARS
{
  valarray<double> f2g(_npoints[dataSetID]);
  F2gamma(dataSetID, f2g, err);
  
  valarray<double> f2gZ(_npoints[dataSetID]);
  F2gammaZ(dataSetID, f2gZ, err);
  
  valarray<double> f2Z(_npoints[dataSetID]);
  F2Z(dataSetID, f2Z, err);      

  valarray<double> k(_npoints[dataSetID]);
  kappa(dataSetID, k);
  // combine together:
  
  double pol     = GetPolarisation(dataSetID);
  double charge = GetCharge(dataSetID);
 
  val = f2g - (_ve + charge*pol*_ae)*k*f2gZ  + (_ae*_ae + _ve*_ve + 2*charge*pol*_ae*_ve)*k * k * f2Z;  
}

void ReactionBaseDISNC::FLgamma BASE_PARS
{
  valarray<double> flu, fld;
  GetFLud(dataSetID, flu,  fld);
  val = 4./9. * flu + 1./9. * fld;
}

void ReactionBaseDISNC::FLgammaZ BASE_PARS
{
  valarray<double> flu, fld;
  GetFLud(dataSetID, flu,  fld);
  val = 2.*( 2./3.*_vu * flu - 1./3.*_vd * fld );
}

void ReactionBaseDISNC::FLZ BASE_PARS
{
  valarray<double> flu, fld;
  GetFLud(dataSetID, flu,  fld);
  val = (_vu*_vu + _au*_au) * flu + (_vd*_vd + _ad*_ad) * fld ;
}


void ReactionBaseDISNC::FL BASE_PARS
{
  valarray<double> flg(_npoints[dataSetID]);
  FLgamma(dataSetID, flg, err);
  
  valarray<double> flgZ(_npoints[dataSetID]);
  FLgammaZ(dataSetID, flgZ, err);
  
  valarray<double> flZ(_npoints[dataSetID]);
  FLZ(dataSetID, flZ, err);      

  valarray<double> k(_npoints[dataSetID]);
  kappa(dataSetID, k);
  // combine together:
  
  double pol     = GetPolarisation(dataSetID);
  double charge = GetCharge(dataSetID);
 
  val = flg - (_ve + charge*pol*_ae)*k*flgZ  + (_ae*_ae + _ve*_ve + 2*charge*pol*_ae*_ve)*k * k * flZ;  
}

void ReactionBaseDISNC::xF3gammaZ BASE_PARS 
{
  valarray<double> xf3u, xf3d;
  GetxF3ud(dataSetID, xf3u,  xf3d);
  val = 2.*(2./3. * _au * xf3u - 1./3. * _ad * xf3d) ;
}

void ReactionBaseDISNC::xF3Z      BASE_PARS 
{
  valarray<double> xf3u, xf3d;
  GetxF3ud(dataSetID, xf3u,  xf3d);
  val = 2.*(_vu * _au * xf3u + _vd * _ad * xf3d) ;
}

void ReactionBaseDISNC::xF3       BASE_PARS 
{
  valarray<double> xf3gZ(_npoints[dataSetID]);
  xF3gammaZ(dataSetID, xf3gZ, err);
  
  valarray<double> xf3Z(_npoints[dataSetID]);
  xF3Z(dataSetID, xf3Z, err);      

  valarray<double> k(_npoints[dataSetID]);
  kappa(dataSetID, k);

  double pol     = GetPolarisation(dataSetID);
  double charge  = GetCharge(dataSetID);

  val =  (_ae*charge + pol*_ve)*k * xf3gZ + (-2*_ae*_ve*charge - pol*(_ve*_ve + _ae*_ae))*k*k * xf3Z;
}

void ReactionBaseDISNC::sred BASE_PARS
{
  auto *yp  = GetBinValues(dataSetID,"y");
  auto y = *yp;

  valarray<double> f2(_npoints[dataSetID]);
  F2(dataSetID, f2, err);

  valarray<double> fl(_npoints[dataSetID]);
  FL(dataSetID, fl, err);

  valarray<double> xf3(_npoints[dataSetID]);
  xF3(dataSetID, xf3, err);

//  double charge = GetCharge(dataSetID);   xF3 is alredy charge-dependent.

  valarray<double> yplus  = 1.0+(1.0-y)*(1.0-y);
  valarray<double> yminus = 1.0-(1.0-y)*(1.0-y);

  val = f2 - y*y/yplus*fl + (yminus/yplus)*xf3 ;  
}


void ReactionBaseDISNC::GetF2ud(int dataSetID, valarray<double>& f2u, valarray<double>& f2d)
{
  // Check if already computed:
  if ( (_f2u[dataSetID])[0] < -99. ) { // compute
  // Get x,Q2 arrays:
    auto *q2p  = GetBinValues(dataSetID,"Q2"), *xp  = GetBinValues(dataSetID,"x");
    auto q2 = *q2p, x = *xp;
    
  // Call QCDNUM
    const int id = 2; const int flag = 0; int Npnt = GetNpoint(dataSetID);
    switch ( GetDataFlav(dataSetID) )
      {
      case dataFlav::incl :
	zmstfun_(id,CNEP2F[0], x[0], q2[0], (_f2u[dataSetID])[0], Npnt, flag);
	zmstfun_(id,CNEM2F[0], x[0], q2[0], (_f2d[dataSetID])[0], Npnt, flag);    
	break ;
      case dataFlav::c :
	zmstfun_(id,CNEP2Fc[0], x[0], q2[0], (_f2u[dataSetID])[0], Npnt, flag);
	break ;
      case dataFlav::b :
	zmstfun_(id,CNEM2Fb[0], x[0], q2[0], (_f2d[dataSetID])[0], Npnt, flag);
	break ;
      }
  }
  f2u = _f2u[dataSetID];
  f2d = _f2d[dataSetID];
}

void ReactionBaseDISNC::GetFLud(int dataSetID, valarray<double>& flu, valarray<double>& fld)
{
  // Check if already computed:
  if ( (_flu[dataSetID])[0] <-99. ) { // compute
    // Get x,Q2 arrays:
    auto *q2p  = GetBinValues(dataSetID,"Q2"), *xp  = GetBinValues(dataSetID,"x");
    auto q2 = *q2p, x = *xp;
    
    // Call QCDNUM
    const int id = 1; const int flag = 0; int Npnt = GetNpoint(dataSetID);
    switch ( GetDataFlav(dataSetID) )
      {
      case dataFlav::incl :
	zmstfun_(id,CNEP2F[0], x[0], q2[0], (_flu[dataSetID])[0], Npnt, flag);
	zmstfun_(id,CNEM2F[0], x[0], q2[0], (_fld[dataSetID])[0], Npnt, flag);    
	break ;
      case dataFlav::c :
	zmstfun_(id,CNEP2Fc[0], x[0], q2[0], (_flu[dataSetID])[0], Npnt, flag);
	break ;     
      case dataFlav::b :
	zmstfun_(id,CNEM2Fb[0], x[0], q2[0], (_fld[dataSetID])[0], Npnt, flag);    
	break ;      
      }
  }
  flu = _flu[dataSetID];
  fld = _fld[dataSetID];

}

void ReactionBaseDISNC::GetxF3ud( int dataSetID, valarray<double>& xf3u, valarray<double>& xf3d )
{
  // Check if already computed:
  if ( (_xf3u[dataSetID])[0] < -99. ) { // compute
    // Get x,Q2 arrays:
    auto *q2p  = GetBinValues(dataSetID,"Q2"), *xp  = GetBinValues(dataSetID,"x");
    auto q2 = *q2p, x = *xp;
    
    // Call QCDNUM
    const int id = 3; const int flag = 0; int Npnt = GetNpoint(dataSetID);
    // OZ 19.10.2017 TODO: F3 is 0 in VFNS for heavy quarks?
    //if ( GetDataType(dataSetID) == dataType::sigred ) {
    if ( GetDataFlav(dataSetID) == dataFlav::incl ) {
      zmstfun_(id,CNEP3F[0], x[0], q2[0], (_xf3u[dataSetID])[0], Npnt, flag);
      zmstfun_(id,CNEM3F[0], x[0], q2[0], (_xf3d[dataSetID])[0], Npnt, flag);    
    }
    else
    {
      for(int p = 0; p < Npnt; p++)
        _xf3u[dataSetID][p] = _xf3d[dataSetID][p] = 0;
    }
    //else {
    //  NOT_IMPLEMENTED("xF3 c,b");
    //}
  }
  xf3u = _xf3u[dataSetID];
  xf3d = _xf3d[dataSetID];
}


void ReactionBaseDISNC::kappa(int dataSetID, valarray<double>& k) {
  auto *q2p = GetBinValues(dataSetID,"Q2");
  double cos2thetaW = 1-_sin2thetaW;

  k= 1./(4*_sin2thetaW*cos2thetaW)  * (*q2p)/( (*q2p)+_Mz*_Mz);  
}
