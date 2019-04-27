
/*
   @file ReactionBaseDISCC.cc
   @date 2017-10-05
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-10-05
*/

#include "ReactionBaseDISCC.h"
#include <iostream>
#include  "QCDNUM/QCDNUM.h"
#include <IntegrateDIS.h>
#include"hf_errlog.h"


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
//TODO: move this to base class
enum class dataFlav{incl,c};        //!< Define final state.
struct ReactionData{
  int          _npoints;                //!< Number of points in a dataset.
  double       _polarisation=0.;        //!< longitudinal polarisation
  double       _charge=0.;              //!< lepton beam charge
  bool         _isReduced=false;        //!< reduced cross section
  dataFlav     _dataFlav=dataFlav::incl;//!< flavour (incl, c, b)
  // for integrated cross sections
  // method is based on legacy subroutine GetIntegratedDisXsection
  IntegrateDIS*_integrated;
  // Some buffering mechanism to avoid double calls
  valarray<double>_f2u; //!< F2 for u-type quarks
  valarray<double>_f2d; //!< F2 for d-type quarks
  valarray<double>_flu; //!< FL for u-type quarks
  valarray<double>_fld; //!< FL for d-type quarks
  valarray<double>_xf3u;
  valarray<double>_xf3d;
};
// Initialize at the start of the computation
void ReactionBaseDISCC::atStart()
{
  // This we do not want to fit:
  _Gf=*XFITTER_PARS::getParamD("gf");
  _convfac=*XFITTER_PARS::getParamD("convFac");
    ///
  int nwords;
  QCDNUM::zmfillw(nwords);//TODO: will this crash if QCDNUM is not initialized?
}

// Main function to compute results at an iteration
void ReactionBaseDISCC::compute(TermData*td,valarray<double>&valExternal,map<string,valarray<double> >&errExternal)
{
  valarray<double> val;//TODO Is this needed?
  map<string, valarray<double> > err;//TODO remove

  // Basic formulae for CC cross section:

 //WIP TermData
  const valarray<double>&y=td->getBinColumn("y");

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

  //for(size_t i = 0; i < f2.size(); i++)
  //  printf("%f %f    %f    %f    %f  =  %f\n", (*GetBinValues(dataSetID,"Q2"))[i], (*GetBinValues(dataSetID,"x"))[i], f2[i], fl[i], xf3[i], val[i]);

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

void ReactionBaseDISCC::atIteration() {
  // Get some basic parameters:
  //huh???
  _MW = GetParam("Mw");

  // Re-set internal maps (faster access):
  // what??? --Ivan
  for ( auto ds : _dsIDs)  {
    (_f2u[ds])[0] = -100.;
    (_flu[ds])[0] = -100.;
    (_xf3u[ds])[0] = -100.;
    (_f2d[ds])[0] = -100.;
    (_fld[ds])[0] = -100.;
    (_xf3d[ds])[0] = -100.;
  }
}
void ReactionBaseDISCC::initTerm(TermData*td){
  ReactionData*rd=new ReactionData();
  td->reactionData=(void*)rd;
  auto&_polarisation=rd->_polarisation;
  auto&_charge      =rd->_charge;
  auto&_isReduced   =rd->_isReduced;
  auto&_integrated  =rd->_integrated;
  if(td->hasParam("epolarity"))_polarisation=*td->getParamD("epolarity");
  if(td->hasParam("echarge"))  _charge      =*td->getParamD("echarge");
  if(td->hasParam("reduced"))  _isReduced   = td->getParamI("reduced");
  // check if settings are provided in the new format key=value
  // type: sigred, signonred (no F2, FL implemented so far, thus type is defined by bool _isReduced)
  // HERA data files provide 'signonred' CC cross sections
  // Inclusive "non-reduced" cross section by default.
  string msg = "I: Calculating DIS CC reduced cross section";
  if(td->hasParam("type")){
    string type=td->getParamS("type");
    if(type=="sigred"){
      _isReduced=true;
      msg = "I: Calculating DIS CC reduced cross section";
    }else if(type=="signonred"){
      _isReduced=false;
      msg = "I: Calculating DIS CC non-reduced cross section";
    }else{
      cerr<<"[ERROR] Unknown type=\""<<type<<"\" given to reaction \""<<getReactionName()<<"\"; termID="<<td->id<<endl;
      hf_errlog(17101903,"F: Unknown \"type\" given to reaction, see stderr");
    }
  }
  // flav: incl, c, b
  if(td->hasParam("flav")){
    string flavor=td->getParamS("flav");
    if(flavor=="incl"){
      _dataFlav=dataFlav::incl;
      msg += " inclusive";
    }else if(flavor=="c"){
      _dataFlav=dataFlav::c;
      msg += " charm";
    }else if(flavor=="b"){//no beauty
      //NOT IMPLEMENTED
      hf_errlog(18042501,"F: predictions for beauty in CC are not available (term id = "+to_string(td->id)+")");
    }else{
      cerr<<"[ERROR] Unknown flavor=\""<<flavor<<"\" given to reaction \""<<getReactionName()<<"\"; termID="<<td->id<<endl;
      hf_errlog(18042502,"F: Unknown \"flavor\" given to reaction, see stderr");
    }
  }
  // check if centre-of-mass energy is provided
  double s = -1.0;
  if(td->hasParam("energy")){
    double energy=*td->getParamD("energy");
    s=energy*energy;
  }
  //TODO: if energy not provided?

  // bins
  // if Q2min, Q2max, ymin and ymax (and optionally xmin, xmax) are provided, integrated cross sections are calculated
  auto*q2minp=td->getBinColumnOrNull("Q2min");
  auto*q2maxp=td->getBinColumnOrNull("Q2max");
  // also try small first letter for Q2 (for backward compatibility)
  if(!q2minp)
       q2minp=td->getBinColumnOrNull("q2min");
  if(!q2maxp)
       q2maxp=td->getBinColumnOrNull("q2max");
  auto*yminp =td->getBinColumnOrNull("ymin");
  auto*ymaxp =td->getBinColumnOrNull("ymax");
  // optional xmin, xmax for integrated cross sections
  auto*xminp =td->getBinColumnOrNull("xmin");
  auto*xmaxp =td->getBinColumnOrNull("xmax");

  if(q2minp && q2maxp && yminp && ymaxp)
  {
    // integrated cross section
    if(s < 0)
      hf_errlog(18060100, "F: centre-of-mass energy is required for integrated DIS dataset " + std::to_string(dataSetID));
    if(_isReduced)
      hf_errlog(18060200, "F: integrated DIS can be calculated only for non-reduced cross sections, dataset " + std::to_string(dataSetID));
    IntegrateDIS* iDIS = new IntegrateDIS();
    _npoints=iDIS->init(s, q2minp, q2maxp, yminp, ymaxp, xminp, xmaxp);
    _integrated=iDIS;
    msg += " (integrated)";
  }
  else
  {
    // cross section at (Q2,x) points
    //TODO: replace with has?
    auto*q2p=&td->getBinColumn("Q2");
    auto*xp =&td->getBinColumn("x");
    auto*yp = td->getBinColumnOrNull("y");

    // if Q2 and x bins and centre-of-mass energy provided, calculate y = Q2 / (s * x)
    if(!yp)
    {
      if ( s > 0.0 )
      {
        valarray<double> y = (*q2p) / (s * (*xp));
        std::pair<string,valarray<double>* > dsBin = std::make_pair("y", &y);
        //Oooof, this is not implemented...
        AddBinning(dataSetID, &dsBin);//TODO: think how to implement this...
        //skipping for now...
        yp = GetBinValues(dataSetID, "y");
      }
      //TODO: else? when s<=0?
    }
    _npoints=q2p->size();
  }

  hf_errlog(17041001,msg);

  // Allocate internal arrays:
  rd->_f2u .resize(_npoints);
  rd->_f2d .resize(_npoints);
  rd->_flu .resize(_npoints);
  rd->_fld .resize(_npoints);
  rd->_xf3u.resize(_npoints);
  rd->_xf3d.resize(_npoints);
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
//TODO: rewrite those
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

