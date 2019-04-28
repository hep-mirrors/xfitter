
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
  //key, x, q2, sf are arrays, better use pointer rather than reference
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
  IntegrateDIS*_integrated=nullptr;
  // Some buffering mechanism to avoid double calls
  valarray<double>_f2u; //!< F2 for u-type quarks
  valarray<double>_f2d; //!< F2 for d-type quarks
  valarray<double>_flu; //!< FL for u-type quarks
  valarray<double>_fld; //!< FL for d-type quarks
  valarray<double>_xf3u;
  valarray<double>_xf3d;
  const double*Mw;//parameter of W mass
  //pointers to bin arrays
  //WIP
  const valarray<double>*q2p;
  const valarray<double>*xp;
  const valarray<double>*yp;
  bool ownYp=false;//true if y-bins were calculated from x and Q2, and *yp was created
};
// Initialize at the start of the computation
void ReactionBaseDISCC::atStart(){
  // This we do not want to fit:
  _Gf     =*XFITTER_PARS::getParamD("gf");
  _convfac=*XFITTER_PARS::getParamD("convFac");
    ///
  int nwords;
  QCDNUM::zmfillw(nwords);//TODO: will this crash if QCDNUM is not initialized?
}

// Main function to compute results at an iteration
void ReactionBaseDISCC::compute(TermData*td,valarray<double>&valExternal,map<string,valarray<double> >&errExternal)
{
  ReactionData*rd=(ReactionData*)td->reactionData;
  const double MW=*rd->Mw;
  size_t _npoints=rd->_npoints;

  // Basic formulae for CC cross section:
  const valarray<double>&y=td->getBinColumn("y");

  valarray<double> yplus  = 1.0+(1.0-y)*(1.0-y);
  valarray<double> yminus = 1.0-(1.0-y)*(1.0-y);

  valarray<double>f2 =F2 (td);
  valarray<double>fl =FL (td);
  valarray<double>xf3=xF3(td);
  double polarity=rd->_polarisation;
  double charge=  rd->_charge;

  valarray<double>val;
  if(charge>0)val=0.5*(1+polarity)*(yplus*f2 - yminus*xf3 - y*y*fl);
  else        val=0.5*(1-polarity)*(yplus*f2 + yminus*xf3 - y*y*fl);

  if(!rd->_isReduced){
    // extra factor for non-reduced cross section
    auto&x =td->getBinColumn("x"),
        &q2=td->getBinColumn("Q2");
    const double pi = 3.1415926535897932384626433832795029;
    valarray<double> factor = (MW*MW*MW*MW/pow((q2+MW*MW),2))*_Gf*_Gf/(2*pi*x)*_convfac;
    val *= factor;
  }

  IntegrateDIS*iDIS=rd->_integrated;
  if(iDIS){
    // integrated cross sections
    valExternal=iDIS->compute(val);
  }else{
    // usual cross section at (q2,x) points
    valExternal=val;
  }
}

  /*
   * I do not understand how that could speed up anything
   * --Ivan
void ReactionBaseDISCC::atIteration() {
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
  */
void ReactionBaseDISCC::initTerm(TermData*td){
  ReactionData*rd=new ReactionData();
  td->reactionData=(void*)rd;
  auto&_isReduced   =rd->_isReduced;
  if(td->hasParam("epolarity"))rd->_polarisation=*td->getParamD("epolarity");//cannot be fitted
  if(td->hasParam("echarge"))  rd->_charge      =*td->getParamD("echarge");  //cannot be fitted
  if(td->hasParam("reduced"))  _isReduced       = td->getParamI("reduced");
  rd->Mw=td->getParamD("Mw");

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
    rd->_integrated=iDIS;
    msg += " (integrated)";
  }
  else
  {
    // cross section at (Q2,x) points
    //TODO: replace with has?
    //Do not use GetBinValues
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
        //maybe I should store pointers in ReactionData and create a new valarray there. That would be best
        //skipping for now...
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

const valarray<double>*GetBinValues(TermData*td,const string&binName){
  IntegrateDIS*iDIS=((ReactionData*)td->reactionData)->_integrated;
  if(iDIS==nullptr)return&td->getBinColumn(binName);
  if     (binName=="Q2")return iDIS->getBinValuesQ2();
  else if(binName=="x" )return iDIS->getBinValuesX();
  else if(binName=="y" )return iDIS->getBinValuesY();
  return&td->getBinColumn(binName);
};

valarray<double>GetF(TermData*td,const int id){
  //other values are illegal
  //return through out
/*F   array id incl   c       isNegative nid
  FL  _flu  1  CCEP2F CCEP2Fc F          0
  F2  _f2u  2  CCEP2F CCEP2Fc F          1
  xF3 _xf3u 3  CCEP3F CCEP3Fc F          2
  FL  _fld  1  CCEM2F CCEM2Fc T          3
  F2  _f2d  2  CCEM2F CCEM2Fc T          4
  xF3 _xf3d 3  CCEM3F CCEM3Fc T          5
*/
  ReactionData*rd=(ReactionData*)td->reactionData;
  bool isNegative=rd->_charge<0;
  char nid=id-1;
  if(isNegative)nid+=3;
  const valarray<double>*f;
  switch(nid){
  case 0:
    f=&rd->_flu;
    break;
  case 1:
    f=&rd->_f2u;
    break;
  case 2:
    f=&rd->_xf3u;
    break;
  case 3:
    f=&rd->_fld;
    break;
  case 4:
    f=&rd->_f2d;
    break;
  case 5:
    f=&rd->_xf3d;
    break;
  default:
    std::abort();//unreachable
  }
  // Check if already computed://I think I broke it. Wait, how did this work, again? --Ivan
  if((*f)[0]<-99){//then compute
    auto&q2=td->getBinColumn("Q2"),
        &x =td->getBinColumn("x");
    const double*C;
    switch(rd->_dataFlav){
    case dataFlav::incl:
      if(isNegative){
        if(id==3)C=CCEM3F;
        else     C=CCEM2F;
      }else{
        if(id==3)C=CCEP3F;
        else     C=CCEP2F;
      }break;
    case dataFlav::c:
      if(isNegative){
        if(id==3)C=CCEM3Fc;
        else     C=CCEM2Fc;
      }else{
        if(id==3)C=CCEP3Fc;
        else     C=CCEP2Fc;
      }break;
    default:
      std::abort();//unreachable
    }
    // Call QCDNUM
    const int flag=0;
    const int Npnt=f->size();
    zmstfun_(id,C[0],x[0],q2[0],(*f)[0],Npnt,flag);
  }
  return*f;
}
valarray<double> FL(TermData*td){return GetF(td,1);}
valarray<double> F2(TermData*td){return GetF(td,2);}
valarray<double>xF3(TermData*td){return GetF(td,3);}
