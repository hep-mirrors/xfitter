 
/*
   @file ReactionBaseDISCC.cc
   @date 2017-10-05
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-10-05
*/

#include "ReactionBaseDISCC.h"
#include <iostream>

// Helpers for QCDNUM (CC): 
const double  CCEP2F[] = {0.,0.,1.,0.,1.,0.,0.,1.,0.,1.,0.,0.,0.} ; 
const double  CCEM2F[] = {0.,0.,0.,1.,0.,1.,0.,0.,1.,0.,1.,0.,0.} ;

const double  CCEP3F[] = {0.,0.,-1.,0.,-1.,0.,0.,1.,0.,1.,0.,0.,0.};
const double  CCEM3F[] = {0.,0. ,0.,-1.,0.,-1.,0.,0.,1.,0.,1.,0.,0.};

// define QCDNUM function:
extern "C" {
  void zmstfun_(const int& id, const double& key, double& x, double& q2, double& sf, const int& np, const int &flag);
}


// the class factories
extern "C" ReactionBaseDISCC* create() {
  return new ReactionBaseDISCC();
}


// Initialize at the start of the computation
int ReactionBaseDISCC::initAtStart(const string &s)
{
  // This we do not want to fit:
  _Gf = GetParam("gf");
  _convfac = GetParam("convFac");
  return 0;
}

// Main function to compute results at an iteration
int ReactionBaseDISCC::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
  // Basic formulat for CC cross section:

  auto *yp  = GetBinValues(dataSetID,"y");
  auto y = *yp;

  valarray<double> yplus  = 1.0+(1.0-y)*(1.0-y);
  valarray<double> yminus = 1.0-(1.0-y)*(1.0-y);


  valarray<double> f2(_npoints[dataSetID]);
  valarray<double> fl(_npoints[dataSetID]);
  valarray<double> xf3(_npoints[dataSetID]);

  F2 (dataSetID,f2,err);
  FL (dataSetID,fl,err);
  xF3(dataSetID,xf3,err);

  double polarity = GetPolarisation(dataSetID);
  
  if ( GetCharge(dataSetID) > 0) {
    val = 0.5*(yplus*f2 - yminus*xf3 - y*y*fl);
    val *= (1+polarity);
  }
  else {
    val = 0.5*(yplus*f2 + yminus*xf3 - y*y*fl);
    val *= (1-polarity);
  }

  if (! IsReduced(dataSetID)) {
    // extra factor for non-reduced cross section
    auto *xp  = GetBinValues(dataSetID,"x");
    auto x = *xp;
    auto *Q2p  = GetBinValues(dataSetID,"Q2");
    auto q2 = *Q2p;
    const double pi = 3.1415926535897932384626433832795029;
    valarray<double> factor = (_MW*_MW*_MW*_MW/pow((q2+_MW*_MW),2))*_Gf*_Gf/(2*pi*x)*_convfac;
    val *= factor;
  }
  
  return 0;
}

void ReactionBaseDISCC::initAtIteration() {
  // Get some basic parameters:
  _MW = GetParam("Mw");
  
  // Re-set internal maps (faster access):
  for ( auto ds : _dsIDs)  {
    (_f2u[ds])[0] = -100.;
    (_flu[ds])[0] = -100.;
    (_xf3u[ds])[0] = -100.;
    (_f2d[ds])[0] = -100.;
    (_fld[ds])[0] = -100.;
    (_xf3d[ds])[0] = -100.;
  }
}

// 
void  ReactionBaseDISCC::setDatasetParamters( int dataSetID, map<string,string> pars, map<string,double> parsDataset) 
{
  auto *q2p  = GetBinValues(dataSetID,"Q2"), *xp  = GetBinValues(dataSetID,"x"), *yp  = GetBinValues(dataSetID,"y");  
  if (q2p == nullptr || xp == nullptr || yp == nullptr ) {
    string msg = "F: Q2, x or Y bins are missing for CC DIS reaction for dataset " + std::to_string(dataSetID);
    hf_errlog_(17100801,msg.c_str(), msg.size());
  }
  _npoints[dataSetID] = (*q2p).size();
  _polarisation[dataSetID] =  (parsDataset.find("epolarity") != parsDataset.end()) ? parsDataset["epolarity"] : 0;
  _charge[dataSetID]       =  (parsDataset.find("echarge")       != parsDataset.end()) ? parsDataset["echarge"] : 0;
  _isReduced[dataSetID]    =  (parsDataset.find("reduced")       != parsDataset.end()) ? parsDataset["reduced"] : 0;

  // check if settings are provided in the new format key=value

  // type: sigred (no F2, FL implemented so far)
  map<string,string>::iterator it = pars.find("type");
  if ( it != pars.end() ) {
    if(it->second == "sigred")
    {
      _isReduced[dataSetID] = 1;
    }
    else
    {
      char buffer[256];
      sprintf(buffer, "F: dataset with id = %d has unknown type = %s", dataSetID, it->second.c_str());
      string str = buffer;
      hf_errlog_(17101903, str.c_str(), str.length());
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

  // Allocate internal arrays:
  _f2u[dataSetID].resize(_npoints[dataSetID]);
  _f2d[dataSetID].resize(_npoints[dataSetID]);
  _flu[dataSetID].resize(_npoints[dataSetID]);
  _fld[dataSetID].resize(_npoints[dataSetID]);
  _xf3u[dataSetID].resize(_npoints[dataSetID]);
  _xf3d[dataSetID].resize(_npoints[dataSetID]);
}


// Get SF

void ReactionBaseDISCC::F2 BASE_PARS
{
  valarray<double> f2;
  if (GetCharge(dataSetID) > 0) {
    GetF2u(dataSetID, f2);
  }
  else {
    GetF2d(dataSetID, f2);
  }
  val = f2;
}

void ReactionBaseDISCC::FL BASE_PARS
{
  valarray<double> fl;
  if (GetCharge(dataSetID) > 0) {
    GetFLu(dataSetID, fl);
  }
  else {
    GetFLd(dataSetID, fl);
  }
  val = fl;
}

void ReactionBaseDISCC::xF3 BASE_PARS
{
  valarray<double> xf3;
  if (GetCharge(dataSetID) > 0) {
    GetxF3u(dataSetID, xf3);
  }
  else {
    GetxF3d(dataSetID, xf3);
  }
  val = xf3;
}


//// -------------------------------------

void ReactionBaseDISCC::GetF2u(int dataSetID, valarray<double>& f2u)
{
  // Check if already computed:
  if ( (_f2u[dataSetID])[0] < -99. ) { // compute
  // Get x,Q2 arrays:
    auto *q2p  = GetBinValues(dataSetID,"Q2"), *xp  = GetBinValues(dataSetID,"x");
    auto q2 = *q2p, x = *xp;
    
  // Call QCDNUM
    const int id = 2; const int flag = 0; int Npnt = GetNpoint(dataSetID);
    zmstfun_(id,CCEP2F[0], x[0], q2[0], (_f2u[dataSetID])[0], Npnt, flag);
    //  zmstfun_(id,CCEM2F[0], x[0], q2[0], (_f2d[dataSetID])[0], Npnt, flag);    
  }
  f2u = _f2u[dataSetID];
  //  f2d = _f2d[dataSetID];
}

void ReactionBaseDISCC::GetFLu(int dataSetID, valarray<double>& flu)
{
  // Check if already computed:
  if ( (_flu[dataSetID])[0] <-99. ) { // compute
    // Get x,Q2 arrays:
    auto *q2p  = GetBinValues(dataSetID,"Q2"), *xp  = GetBinValues(dataSetID,"x");
    auto q2 = *q2p, x = *xp;
    
    // Call QCDNUM
    const int id = 1; const int flag = 0; int Npnt = GetNpoint(dataSetID);
    zmstfun_(id,CCEP2F[0], x[0], q2[0], (_flu[dataSetID])[0], Npnt, flag);
    //    zmstfun_(id,CCEM2F[0], x[0], q2[0], (_fld[dataSetID])[0], Npnt, flag);    
  }
  flu = _flu[dataSetID];
  // fld = _fld[dataSetID];

}

void ReactionBaseDISCC::GetxF3u( int dataSetID, valarray<double>& xf3u )
{
  // Check if already computed:
  if ( (_xf3u[dataSetID])[0] < -99. ) { // compute
    // Get x,Q2 arrays:
    auto *q2p  = GetBinValues(dataSetID,"Q2"), *xp  = GetBinValues(dataSetID,"x");
    auto q2 = *q2p, x = *xp;
    
    // Call QCDNUM
    const int id = 3; const int flag = 0; int Npnt = GetNpoint(dataSetID);

    zmstfun_(id,CCEP3F[0], x[0], q2[0], (_xf3u[dataSetID])[0], Npnt, flag);
    //    zmstfun_(id,CCEM3F[0], x[0], q2[0], (_xf3d[dataSetID])[0], Npnt, flag);    
  }
  xf3u = _xf3u[dataSetID];
  //xf3d = _xf3d[dataSetID];
}


void ReactionBaseDISCC::GetF2d(int dataSetID, valarray<double>& f2d)
{
  // Check if already computed:
  if ( (_f2d[dataSetID])[0] < -99. ) { // compute
  // Get x,Q2 arrays:
    auto *q2p  = GetBinValues(dataSetID,"Q2"), *xp  = GetBinValues(dataSetID,"x");
    auto q2 = *q2p, x = *xp;
    
  // Call QCDNUM
    const int id = 2; const int flag = 0; int Npnt = GetNpoint(dataSetID);
    // zmstfun_(id,CCEP2F[0], x[0], q2[0], (_f2u[dataSetID])[0], Npnt, flag);
    zmstfun_(id,CCEM2F[0], x[0], q2[0], (_f2d[dataSetID])[0], Npnt, flag);    
  }
  //f2u = _f2u[dataSetID];
  f2d = _f2d[dataSetID];
}

void ReactionBaseDISCC::GetFLd(int dataSetID, valarray<double>& fld)
{
  // Check if already computed:
  if ( (_fld[dataSetID])[0] <-99. ) { // compute
    // Get x,Q2 arrays:
    auto *q2p  = GetBinValues(dataSetID,"Q2"), *xp  = GetBinValues(dataSetID,"x");
    auto q2 = *q2p, x = *xp;
    
    // Call QCDNUM
    const int id = 1; const int flag = 0; int Npnt = GetNpoint(dataSetID);
    // zmstfun_(id,CCEP2F[0], x[0], q2[0], (_flu[dataSetID])[0], Npnt, flag);
    zmstfun_(id,CCEM2F[0], x[0], q2[0], (_fld[dataSetID])[0], Npnt, flag);    
  }
  // flu = _flu[dataSetID];
  fld = _fld[dataSetID];

}

void ReactionBaseDISCC::GetxF3d( int dataSetID, valarray<double>& xf3d )
{
  // Check if already computed:
  if ( (_xf3d[dataSetID])[0] < -99. ) { // compute
    // Get x,Q2 arrays:
    auto *q2p  = GetBinValues(dataSetID,"Q2"), *xp  = GetBinValues(dataSetID,"x");
    auto q2 = *q2p, x = *xp;
    
    // Call QCDNUM
    const int id = 3; const int flag = 0; int Npnt = GetNpoint(dataSetID);

    // zmstfun_(id,CCEP3F[0], x[0], q2[0], (_xf3u[dataSetID])[0], Npnt, flag);
    zmstfun_(id,CCEM3F[0], x[0], q2[0], (_xf3d[dataSetID])[0], Npnt, flag);    
  }
  //xf3u = _xf3u[dataSetID];
  xf3d = _xf3d[dataSetID];
}

