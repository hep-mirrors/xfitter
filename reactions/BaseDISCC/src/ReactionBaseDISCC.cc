 
/*
   @file ReactionBaseDISCC.cc
   @date 2017-10-05
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-10-05
*/

#include "ReactionBaseDISCC.h"
#include <iostream>
#include <IntegrateDIS.h>

// Helpers for QCDNUM (CC): 

//! full
const double  CCEP2F[] = {0.,0.,1.,0.,1.,0., 0. ,1.,0.,1.,0.,0.,0.} ;
const double  CCEM2F[] = {0.,0.,0.,1.,0.,1., 0. ,0.,1.,0.,1.,0.,0.} ;

const double  CCEP3F[] = {0.,0.,-1.,0.,-1.,0., 0., 1.,0.,1.,0.,0.,0.};
const double  CCEM3F[] = {0.,0. ,0.,-1.,0.,-1., 0., 0.,1.,0.,1.,0.,0.};

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
const double  CCEP2Fc[] = {0.,0.,1.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.} ;
const double  CCEM2Fc[] = {0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,1.,0.,0.} ;

// only c
//const double  CCEP3Fc[] = {0.,0.,-1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
//const double  CCEM3Fc[] = {0.,0. ,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.};
// only s
//const double  CCEP3Fc[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.};
//const double  CCEM3Fc[] = {0.,0.,0.,-1.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
// only s,c
const double  CCEP3Fc[] = {0.,0.,-1.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.};
const double  CCEM3Fc[] = {0.,0.,0.,-1.,0.,0.,0.,0.,0.,0.,1.,0.,0.};

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
int ReactionBaseDISCC::compute(int dataSetID, valarray<double> &valExternal, map<string, valarray<double> > &errExternal)
{
  valarray<double> val;
  map<string, valarray<double> > err;

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
    if(_stFun[dataSetID] == stFun::all)
      val = 0.5*(yplus*f2 - yminus*xf3 - y*y*fl);
    else if(_stFun[dataSetID] == stFun::f2)
    {
      if(IsReduced(dataSetID))
        val = 0.5*(yplus*f2);
      else
        val = f2;
    }
    else if(_stFun[dataSetID] == stFun::fl)
    {
      if(IsReduced(dataSetID))
        val = 0.5*( - y*y*fl);
      else
        val = fl;
    }
    else if(_stFun[dataSetID] == stFun::xf3)
    {
      if(IsReduced(dataSetID))
        val = 0.5*( - yminus*xf3);
      else
        val = xf3;
    }
    val *= (1+polarity);
  }
  else {
    if(_stFun[dataSetID] == stFun::all)
      val = 0.5*(yplus*f2 + yminus*xf3 - y*y*fl);
    else if(_stFun[dataSetID] == stFun::f2)
    {
      if(IsReduced(dataSetID))
        val = 0.5*(yplus*f2);
      else
        val = f2;
    }
    else if(_stFun[dataSetID] == stFun::fl)
    {
      if(IsReduced(dataSetID))
        val = 0.5*( - y*y*fl);
      else
        val = fl;
    }
    else if(_stFun[dataSetID] == stFun::xf3)
    {
      if(IsReduced(dataSetID))
        val = 0.5*(yminus*xf3);
      else
        val = xf3;
    }
    val *= (1-polarity);
  }

  if(!IsReduced(dataSetID) && _stFun[dataSetID] == stFun::all)
  {
    // transform reduced -> non-reduced (double-differential) cross sections
    auto *xp  = GetBinValues(dataSetID,"x");
    auto x = *xp;
    auto *Q2p  = GetBinValues(dataSetID,"Q2");
    auto q2 = *Q2p;
    const double pi = 3.1415926535897932384626433832795029;
    valarray<double> factor = (_MW*_MW*_MW*_MW/pow((q2+_MW*_MW),2))*_Gf*_Gf/(2*pi*x)*_convfac;
    val *= factor;
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
void  ReactionBaseDISCC::setDatasetParameters( int dataSetID, map<string,string> pars, map<string,double> parsDataset) 
{
  _polarisation[dataSetID] =  (parsDataset.find("epolarity") != parsDataset.end()) ? parsDataset["epolarity"] : 0;
  _charge[dataSetID]       =  (parsDataset.find("echarge")       != parsDataset.end()) ? parsDataset["echarge"] : 0;
  _isReduced[dataSetID]    =  (parsDataset.find("reduced")       != parsDataset.end()) ? parsDataset["reduced"] : 0;

  // check if settings are provided in the new format key=value
  // type: signonred, sigred
  // stfun: all (default), f2, fl, xf3
  // NOTE: for sigred = 0, if stfun = all then double-diff. cross sections is calculated,
  // but if stfun = f2, fl and xf3 then structure function is calculated
  // example: sigred = 1, stfun = all -> calculate reduced cross section
  // example: sigred = 0, stfun = all -> calculate double-differential (non-reduced) cross section
  // example: sigred = 0, stfun = f3 -> calculate xF3
  // example: sigred = 1, stfun = f3 -> calculate xF3 contribution to reduced cross section, i.e. (yminus/yplus)*xf3
  _dataFlav[dataSetID] = dataFlav::incl;
  string msg = "I: Calculating DIS CC reduced cross section";
  map<string,string>::iterator it = pars.find("type");
  if ( it != pars.end() ) {
    if(it->second == "sigred")
    {
      _isReduced[dataSetID] = 1;
      msg = "I: Calculating DIS CC reduced cross section";
    }
    else if(it->second == "signonred")
    {
      _isReduced[dataSetID] = 0;
      msg = "I: Calculating DIS CC non-reduced cross section";
    }
    else
    {
      char buffer[256];
      sprintf(buffer, "F: dataset with id = %d has unknown type = %s", dataSetID, it->second.c_str());
      string str = buffer;
      hf_errlog_(17101903, str.c_str(), str.length());
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
    // no beauty
    else if(it->second == "b")
    {
      char buffer[256];
      sprintf(buffer, "F: predictions for beauty in CC are not available (dataset id = %d)", dataSetID);
      string str = buffer;
      hf_errlog_(18042501, str.c_str(), str.length());
    }
    else
    {
      char buffer[256];
      sprintf(buffer, "F: dataset with id = %d has unknown flav = %s", dataSetID, it->second.c_str());
      string str = buffer;
      hf_errlog_(18042502, str.c_str(), str.length());
    }
  }

  // structrure function contrbution: all (default), f2, fl, f3
  _stFun[dataSetID] = stFun::all;
  it = pars.find("stfun");
  if ( it != pars.end() )
  {
    if(it->second == "f2")
    {
      _stFun[dataSetID] = stFun::f2;
      msg += " (F2)";
    }
    else if(it->second == "fl")
    {
      _stFun[dataSetID] = stFun::fl;
      msg += " (FL)";
    }
    else if(it->second == "xf3")
    {
      _stFun[dataSetID] = stFun::xf3;
      msg += " (F3)";
    }
    else if(it->second == "all")
    {
      // do notinng: default option
    }
    else
    {
      char buffer[256];
      sprintf(buffer, "F: dataset with id = %d has unknown stfun = %s", dataSetID, it->second.c_str());
      string str = buffer;
      hf_errlog_(18042502, str.c_str(), str.length());
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
  // if Q2min, Q2max, ymin and ymax (and optionally xmin, xmax) are provided, integrated cross sections are calculated
  auto *q2minp  = GetBinValues(dataSetID,"Q2min");
  auto *q2maxp  = GetBinValues(dataSetID,"Q2max");
  // also try small first letter for Q2 (for backward compatibility)
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
    if(IsReduced(dataSetID))
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
      string msg = "F: Q2, x or Y bins are missing for CC DIS reaction for dataset " + std::to_string(dataSetID);
      hf_errlog_(17100801,msg.c_str(), msg.size());
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

valarray<double> *ReactionBaseDISCC::GetBinValues(int idDS, const string& binName)
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
    switch ( GetDataFlav(dataSetID) )
      {
      case dataFlav::incl :
        zmstfun_(id,CCEP2F[0], x[0], q2[0], (_f2u[dataSetID])[0], Npnt, flag);
        break;
      case dataFlav::c :
        zmstfun_(id,CCEP2Fc[0], x[0], q2[0], (_f2u[dataSetID])[0], Npnt, flag);
        break ;
      }
    //for(int i = 0; i < Npnt; i++)
    //  printf("%f %f    %f\n", q2[i], x[i], (_f2u[dataSetID])[i]);
  }
  f2u = _f2u[dataSetID];
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
    switch ( GetDataFlav(dataSetID) )
      {
      case dataFlav::incl :
        zmstfun_(id,CCEP2F[0], x[0], q2[0], (_flu[dataSetID])[0], Npnt, flag);
        break;
      case dataFlav::c :
        zmstfun_(id,CCEP2Fc[0], x[0], q2[0], (_flu[dataSetID])[0], Npnt, flag);
        break ;
      }
  }
  flu = _flu[dataSetID];

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
    switch ( GetDataFlav(dataSetID) )
      {
      case dataFlav::incl :
        zmstfun_(id,CCEP3F[0], x[0], q2[0], (_xf3u[dataSetID])[0], Npnt, flag);
        break;
      case dataFlav::c :
        //printf("before XF3u: %f\n", (_xf3u[dataSetID])[_xf3u[dataSetID].size() - 1]);
        zmstfun_(id,CCEP3Fc[0], x[0], q2[0], (_xf3u[dataSetID])[0], Npnt, flag);
        //printf("after XF3u: %f\n", (_xf3u[dataSetID])[_xf3u[dataSetID].size() - 1]);
        break;
      }
    //_xf3u[dataSetID] = _xf3u[dataSetID] * x;
  }
  //printf("XF3u: %f\n", (_xf3u[dataSetID])[_xf3u[dataSetID].size() - 1]);
  xf3u = _xf3u[dataSetID];
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
    switch ( GetDataFlav(dataSetID) )
      {
      case dataFlav::incl :
        zmstfun_(id,CCEM2F[0], x[0], q2[0], (_f2d[dataSetID])[0], Npnt, flag);
        break;
      case dataFlav::c :
        zmstfun_(id,CCEM2Fc[0], x[0], q2[0], (_f2d[dataSetID])[0], Npnt, flag);
        break ;
      }
  }
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
    switch ( GetDataFlav(dataSetID) )
      {
      case dataFlav::incl :
        zmstfun_(id,CCEM2F[0], x[0], q2[0], (_fld[dataSetID])[0], Npnt, flag);
        break;
      case dataFlav::c :
        zmstfun_(id,CCEM2Fc[0], x[0], q2[0], (_fld[dataSetID])[0], Npnt, flag);
        break ;
      }
  }
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
    switch ( GetDataFlav(dataSetID) )
      {
      case dataFlav::incl :
        zmstfun_(id,CCEM3F[0], x[0], q2[0], (_xf3d[dataSetID])[0], Npnt, flag);
        break;
      case dataFlav::c :
        zmstfun_(id,CCEM3Fc[0], x[0], q2[0], (_xf3d[dataSetID])[0], Npnt, flag);
        break;
    }
    //_xf3d[dataSetID] = _xf3d[dataSetID] * x;
  }
  xf3d = _xf3d[dataSetID];
}

