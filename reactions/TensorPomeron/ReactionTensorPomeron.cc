 
/*
   @file ReactionTensorPomeron.cc
   @date 2017-04-14
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-04-14
*/

#include "ReactionTensorPomeron.h"
#include <math.h>
#include <fstream>
#include <string.h>
#include <iostream>
#include "xfitter_pars.h"
#include "hf_errlog.h"

// the class factories
extern "C" ReactionTensorPomeron* create() {
  return new ReactionTensorPomeron();
}

vector<double> GetParamV(const string& name)                                                                                                       
{                                                                                                                                                             
  return XFITTER_PARS::gParametersV.at(name);                                                                                                                          
}                                    

vector<double> ReactionTensorPomeron::getSplinePar(const std::string& vn){
  size_t len =  GetParamV(vn).size();
  vector<double> out(len);
  // fudge the name
  std::string n2 = vn.substr(0,3);

  for (size_t i=0; i<len; i++) {
    out[i] = *XFITTER_PARS::getParamD(n2+"v"+std::to_string(i));
  }
  return out;
}

void ReactionTensorPomeron::atIteration() 
{
  // Get spline pars  

  vector<double> s0bv = getSplinePar("s0bn");
  _s0b.set_points(_s0bn, s0bv);

  vector<double> s1bv = getSplinePar("s1bn");
  _s1b.set_points(_s1bn, s1bv);

  if ( _splineR ) {
    vector<double> s0rv = getSplinePar("s0rn");
    _s0r.set_points(_s0rn, s0rv);

    vector<double> s1rv = getSplinePar("s1rn");
    _s1r.set_points(_s1rn, s1rv);
  }
  else {
    // Different param. for R:
    _r0 = *XFITTER_PARS::getParamD("r0");
    _r1 = *XFITTER_PARS::getParamD("r1");
    _delta0 = *XFITTER_PARS::getParamD("delta0");
    _delta1 = *XFITTER_PARS::getParamD("delta1");
  }

  _epsilon0 = *XFITTER_PARS::getParamD("epsilon0");
  _epsilon1 = *XFITTER_PARS::getParamD("epsilon1");


  // Also reggeon
  _c0ir = *XFITTER_PARS::getParamD("c0ir");
  _c1ir = *XFITTER_PARS::getParamD("c1ir");
  _alphaIR = *XFITTER_PARS::getParamD("alphaIR");
  _epsilonR = _alphaIR - 1.;

  // Update internal valarrays too
  for (int ds : _dsIDs) {
    auto q2 = _q2[ds];
    for ( size_t i=0; i<(q2).size(); i++) {
      _b0[ds][i] = b0q2((q2)[i]);
      _b1[ds][i] = b1q2((q2)[i]);
      _qa0[ds][i] = q2a0((q2)[i]);
      _qa1[ds][i] = q2a1((q2)[i]);

      _b2[ds][i]  = b2q2((q2)[i]);
      _qa2[ds][i] = q2a2((q2)[i]);

    }
  }
}

void  ReactionTensorPomeron::initTerm(TermData* td)
{
  int dataSetID = td->id;
  _dsIDs.push_back(dataSetID);
  
  // First see if ep or gamma p:
  _isReduced[dataSetID] = true;
  if (td->hasParam("type")) {
    string type = td->getParamS("type");
    if (type == "signonred")  {
      _isReduced[dataSetID] = false;
    }
    else if (type == "sigred") {
      _isReduced[dataSetID] = true;
    }
  }
  

  // Proton mass
  _mp   = *td->getParamD("mP");


  if (_isReduced[dataSetID]) {
    auto *q2p  = GetBinValues(td,"Q2"), *xp  = GetBinValues(td,"x"), *yp  = GetBinValues(td,"y");  
    if (q2p == nullptr || xp == nullptr || yp == nullptr ) {      
      hf_errlog(17040801,"F: Q2, x or Y bins are missing for NC DIS reaction for dataset " + std::to_string(dataSetID));
    }
    _npoints[dataSetID] = (*q2p).size();
    // Basic kinematic vars:
    _x[dataSetID] = *xp;
    _q2[dataSetID] = *q2p;
    _y[dataSetID] = *yp;
    _W2[dataSetID] = (*q2p)*(1./(*xp)-1.) + _mp*_mp;
  }
  else {
    auto *wp = GetBinValues(td,"W") ;
    if (wp == nullptr) {
      hf_errlog(17042201,"F: W bins are missing for gamma-p reaction for dataset " + std::to_string(dataSetID));
    }
    _npoints[dataSetID] = (*wp).size();
    _W2[dataSetID] = (*wp)*(*wp);
    _q2[dataSetID].resize(_npoints[dataSetID]);
    for ( auto  &q2 : _q2[dataSetID] ) {
      q2 = 0.0;
    }
    _y[dataSetID].resize(_npoints[dataSetID]);
    _x[dataSetID].resize(_npoints[dataSetID]);

  }

  _s0bn = GetParamV("s0bn");
  _s1bn = GetParamV("s1bn");
  _s0rn = GetParamV("s0rn");
  _s1rn = GetParamV("s1rn");

  _m02 = *td->getParamD("m02");
  _m12 = *td->getParamD("m12");

  _convFact = *td->getParamD("ubarnFromGeV");

  // goto log, add offset:
  auto toLog = []( vector<double>& a, double b) { 
    for (size_t i=0; i<a.size(); i++) {
      a[i] = log(a[i]+b);
    } 
  };

  toLog(_s0bn,_m02); toLog(_s1bn,_m12); toLog(_s0rn,0.0); toLog(_s1rn,0.0);

  _alpha_em = *td->getParamD("alphaEM");

  _alphaP = *td->getParamD("alphaP");
  _beta   = *td->getParamD("beta");

  _splineR = td->getParamI("SplineForR"); 

  if (_splineR) {
    hf_errlog(17042201,"I: Unsing spline param. for R in TensorPomeron model");
  }
  _mr = *td->getParamD("mr");

  _alphaR = *td->getParamD("alphaR");
  _betaR  = *td->getParamD("betaR");

  _pq[dataSetID] = 0.5*(_W2[dataSetID]+_q2[dataSetID]-_mp*_mp);
  _delta[dataSetID] = 0.5*_mp*_mp*_q2[dataSetID]/(_pq[dataSetID]*_pq[dataSetID]);
  _Yp[dataSetID] = 1. + (1. - _y[dataSetID])*(1. - _y[dataSetID]) + _y[dataSetID]*_y[dataSetID]*_delta[dataSetID];
  _f[dataSetID] = _y[dataSetID]*_y[dataSetID]*(2.0*_delta[dataSetID]+1.0)/_Yp[dataSetID];

  _b0[dataSetID].resize(_npoints[dataSetID]);
  _b1[dataSetID].resize(_npoints[dataSetID]);
  _b2[dataSetID].resize(_npoints[dataSetID]);
  _qa0[dataSetID].resize(_npoints[dataSetID]);
  _qa1[dataSetID].resize(_npoints[dataSetID]);
  _qa2[dataSetID].resize(_npoints[dataSetID]);
}

// Main function to compute results at an iteration
void ReactionTensorPomeron::compute(TermData* td, valarray<double> &val, map<string, valarray<double> > &err)
{
  int dataSetID = td->id;
  int Np = _npoints[dataSetID];

  valarray<double> sL(Np);
  valarray<double> sLT(Np);

  auto W2 = _W2[dataSetID] ;
  auto Q2 = _q2[dataSetID] ;
  auto pq = _pq[dataSetID] ;
  auto Yp = _Yp[dataSetID] ;
  auto  y = _y[dataSetID] ;
  auto  f = _f[dataSetID] ;
  auto delta = _delta[dataSetID] ;

  if ( _isReduced[dataSetID] ) {
    sigma_L(dataSetID, sL);
    sigma_LT(dataSetID, sLT);
  
    double prefix           = 1.0/(4.*M_PI*M_PI*_alpha_em);
    valarray<double> yf     = Yp/(1.0+(1.0-y)*(1.0-y));
    valarray<double> df     = 1.0/(1.0+2.0*delta);
    valarray<double> kf     = 0.5 * Q2*(W2 - _mp*_mp)/pq;

    val = prefix*yf*df*kf*(sLT - f*sL);
  }
  else {
    sigma_LT(dataSetID, sLT);
    val = _convFact*sLT;
  }
  return;
}

void ReactionTensorPomeron::sigma_L(int dataSetID, valarray<double>& sL) 
{
  auto W2 = _W2[dataSetID] ;
  auto Q2 = _q2[dataSetID] ;
  auto pq = _pq[dataSetID] ;

  auto b0 = _b0[dataSetID] ;
  auto b1 = _b1[dataSetID] ;
  auto qa0 = _qa0[dataSetID] ;
  auto qa1 = _qa1[dataSetID] ;

  // Reggeon:
  auto b2 = _b2[dataSetID];
  auto qa2 = _qa2[dataSetID];


  valarray<double> prefix =  4. * M_PI * _alpha_em / ( W2*(W2 - _mp*_mp) ) ;
  valarray<double> suf_a =  (4.* pq *pq + Q2*_mp*_mp) ;
  valarray<double> suf_b =  ( Q2*_mp*_mp );

  valarray<double> p0 = 3*_beta * cos(_epsilon0 * M_PI / 2.0) * pow(_alphaP*W2, _epsilon0) ;
  valarray<double> p1 = 3*_beta * cos(_epsilon1 * M_PI / 2.0) * pow(_alphaP*W2, _epsilon1) ;
  
  // Reggeon
  valarray<double> r  = 3*_beta * cos(_epsilonR * M_PI / 2.0) * pow(_alphaP*W2, _epsilonR) ;

  sL = prefix*( 
	       p0*2.0*(qa0*suf_a + b0*suf_b) 
	       + p1*2.0*(qa1*suf_a + b1*suf_b) 
	       + r*2.0*(qa2*suf_a + b2*suf_b) 
		);
}

void ReactionTensorPomeron::sigma_LT(int dataSetID, valarray<double>& sLT) 
{
  auto W2 = _W2[dataSetID] ;
  auto Q2 = _q2[dataSetID] ;
  auto pq = _pq[dataSetID] ;

  auto b0 = _b0[dataSetID] ;
  auto b1 = _b1[dataSetID] ;
  auto b2 = _b2[dataSetID] ;

  valarray<double> prefix =  4. * M_PI * _alpha_em / ( W2*(W2 - _mp*_mp) ) ;

  valarray<double> suffix =  ( pq *pq + Q2*_mp*_mp) ;

  valarray<double> p0 = 3*_beta * cos(_epsilon0 * M_PI / 2.0) * pow(_alphaP*W2, _epsilon0) * 4.* b0;
  valarray<double> p1 = 3*_beta * cos(_epsilon1 * M_PI / 2.0) * pow(_alphaP*W2, _epsilon1) * 4.* b1;
  // Reggeon:
  valarray<double>  r = 3*_beta * cos(_epsilonR * M_PI / 2.0) * pow(_alphaP*W2, _epsilonR) * 4.* b2; 

  sLT = prefix*(p0+p1 +  r )*suffix;
}

void ReactionTensorPomeron::atFCN3() {
  writeOut("pomeron.csv");
  // also write out the parameters, to make next iteration easier
  //  std::cout <<  emitReactionLocalPars() << std::endl;
}

void ReactionTensorPomeron::atMakeErrorBands(int ivector) {
  std::string FileName = "pom_"+std::to_string(ivector)+".csv";
  writeOut(FileName);
}

void ReactionTensorPomeron::writeOut(const std::string& file) {
  std::ofstream f;
  f.open(file) ; // XXXX should go to ./output !!!
  f << "logQ2  q2a0  q2a1  b0  b1  r0  r1 " << std::endl;
  for ( double q2l =-3; q2l<log(50.0); q2l += 0.1) {
    double q2 = exp(q2l);
    f << q2l 
      << " " << q2a0(q2) 
      << " " << q2a1(q2) 
      << " " << b0q2(q2) 
      << " " << b1q2(q2)  
      << " " << r0q2(q2) 
      << " " << r1q2(q2)  
      << std::endl;
  }
  f.close();

}

const double ReactionTensorPomeron::r0q2(double q2) { 
  if (_splineR) {
    return exp( _s0r(log(q2)) ); 
  }
  else {
    return _r0*q2 * pow( (1.0 + q2/(_mr*_mr)), _delta0-1); 
  }
}

const double ReactionTensorPomeron::r1q2(double q2) { 
  if (_splineR) {
    return exp( _s1r(log(q2)) ); 
  }
  else {
    return _r1*q2 * pow( (1.0 + q2/(_mr*_mr)), _delta1-1); 
  }
}

const valarray<double> *ReactionTensorPomeron::GetBinValues(TermData *td, const string &binName)                                                                    
{                                                                                                                                                               
  unsigned termID = td->id;                                                                                                                                     
                                                                                                                                                                
  return td->getBinColumnOrNull(binName);                                                                                                                     
}
