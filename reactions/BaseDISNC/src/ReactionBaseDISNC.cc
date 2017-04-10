 
/*
   @file ReactionBaseDISNC.cc
   @date 2017-04-08
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-04-08
*/

#include "ReactionBaseDISNC.h"
#include <iostream>

template <typename T>
void print(T d) {
  std::cout << d << "\n";
}

// Helpers for QCDNUM:

//! F2,FL full
const double CNEP2F[] = {0.,0.,1.,0.,1.,0.,0.,0.,1.,0.,1.,0.,0.}; //u  (top off ?)
const double CNEM2F[] = {0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0.}; //d

//! xF3 full
const double CNEP3F[] = {0.,0.,-1.,0.,-1.,0.,0.,0.,1.,0.,1.,0.,0.}; //u
const double CNEM3F[] = {0.,-1.,0.,-1.,0.,-1.,0.,1.,0.,1.,0.,1.,0.}; //d

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
  return 0;
}

// Main function to compute results at an iteration
int ReactionBaseDISNC::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
  switch ( _dataType) 
    {
    case dataType::sigred :
      sred(dataSetID, val, err) ;
      break ;
    case dataType::f2c :
      F2(dataSetID, val, err) ;
      break ;
    case dataType::f2b :
      F2(dataSetID, val, err) ;
      break ;
    }
  return 0;
}

void ReactionBaseDISNC::initAtIteration() {
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
void  ReactionBaseDISNC::setDatasetParamters( int dataSetID, map<string,string> pars, map<string,double> parsDataset) 
{
  auto *q2p  = GetBinValues(dataSetID,"Q2"), *xp  = GetBinValues(dataSetID,"x"), *yp  = GetBinValues(dataSetID,"y");  
  if (q2p == nullptr || xp == nullptr || yp == nullptr ) {
    string msg = "F: Q2, x or Y bins are missing for NC DIS reaction for dataset " + std::to_string(dataSetID);
    hf_errlog_(17040801,msg.c_str(), msg.size());
  }
  _npoints[dataSetID] = (*q2p).size();
  _polarisation[dataSetID] =  (parsDataset.find("epolarity") != parsDataset.end()) ? parsDataset["epolarity"] : 0;
  _charge[dataSetID]       =  (parsDataset.find("echarge")       != parsDataset.end()) ? parsDataset["echarge"] : 0;

  _dataType = dataType::sigred;  // Reduced cross section by default.
  string msg = "I: Calculating DIS NC reduced cross section";
  if ( parsDataset.find("F2c") != parsDataset.end() ) {
    _dataType = dataType::f2c;    
    msg = "I: Calculating DIS NC F2c";
  }
  if ( parsDataset.find("F2b") != parsDataset.end() ) {
    _dataType = dataType::f2b;
    msg = "I: Calculating DIS NC F2b";
  }
  if ( parsDataset.find("reduced") != parsDataset.end() ) {
    _dataType = dataType::sigred;
    msg = "I: Calculating DIS NC reduced cross section";
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
  val = 2*( 2./3.*_vu * f2u - 1./3.*_vd * f2d );
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
  val = 2*( 2./3.*_vu * flu - 1./3.*_vd * fld );
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
  val = 2/3 * _au * xf3u - 1/3 * _ad * xf3d ;
}

void ReactionBaseDISNC::xF3Z      BASE_PARS 
{
  valarray<double> xf3u, xf3d;
  GetxF3ud(dataSetID, xf3u,  xf3d);
  val = _vu * _au * xf3u + _vd * _ad * xf3d ;
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

  val = -(_ae + charge*pol*_ve)*k * xf3gZ + (2*_ae*_ve + charge*pol*(_ve*_ve + _ae*_ae))*k*k * xf3Z;
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

  double charge = GetCharge(dataSetID);

  valarray<double> yplus  = 1.0+(1.0-y)*(1.0-y);
  valarray<double> yminus = 1.0-(1.0-y)*(1.0-y);

  val = f2 - y*y/yplus*fl - charge*(yminus/yplus)*xf3 ;

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
    switch ( _dataType )
      {
      case dataType::sigred :
	zmstfun_(id,CNEP2F[0], x[0], q2[0], (_f2u[dataSetID])[0], Npnt, flag);
	zmstfun_(id,CNEM2F[0], x[0], q2[0], (_f2d[dataSetID])[0], Npnt, flag);    
	break ;
      case dataType::f2c :
	zmstfun_(id,CNEP2Fc[0], x[0], q2[0], (_f2u[dataSetID])[0], Npnt, flag);
	break ;
      case dataType::f2b :
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
    switch ( _dataType )
      {
      case dataType::sigred : 
	zmstfun_(id,CNEP2F[0], x[0], q2[0], (_flu[dataSetID])[0], Npnt, flag);
	zmstfun_(id,CNEM2F[0], x[0], q2[0], (_fld[dataSetID])[0], Npnt, flag);    
	break ;
      case dataType::f2c : 
	zmstfun_(id,CNEP2Fc[0], x[0], q2[0], (_flu[dataSetID])[0], Npnt, flag);
	break ;     
      case dataType::f2b : 
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
    if ( _dataType == dataType::sigred ) {
      zmstfun_(id,CNEP3F[0], x[0], q2[0], (_xf3u[dataSetID])[0], Npnt, flag);
      zmstfun_(id,CNEM3F[0], x[0], q2[0], (_xf3d[dataSetID])[0], Npnt, flag);    
    }
    else {
      NOT_IMPLEMENTED(" xF3 b,c ");
    }
  }
  xf3u = _xf3u[dataSetID];
  xf3d = _xf3d[dataSetID];
}


void ReactionBaseDISNC::kappa(int dataSetID, valarray<double>& k) {
  auto *q2p = GetBinValues(dataSetID,"Q2");
  double cos2thetaW = sqrt(1-_sin2thetaW*_sin2thetaW);

  k= 1./(4*_sin2thetaW*cos2thetaW)  * (*q2p)/( (*q2p)+_Mz*_Mz);

}
