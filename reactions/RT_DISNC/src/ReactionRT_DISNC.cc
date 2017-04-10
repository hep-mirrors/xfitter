/*
   @file ReactionRT_DISNC.cc
   @date 2017-04-10
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-04-10
*/

#include <iostream>
#include "ReactionRT_DISNC.h"

// the class factories
extern "C" ReactionRT_DISNC* create() {
  return new ReactionRT_DISNC();
}

// RT wrappers from RT/src/mstw2008_wrap.f:
extern "C" {
  void mstwnc_wrap_(const double& x, const double& q2, const int& ipn, 
		    double& f2, double& f2c, double& f2b,  double& fl,  double& flc, double& flb, 
		    const int& iflag, const int& index, const double& f2QCDNUM, const double& flQCDNUM,
		    const int& usekfactors = 0);
  void rt_setalphas_(const double& alphaSzero);
  void rt_set_input_(const double* varin, const double& mCharmin, const double& mBottomin, const double& alphaSQ0in,
		     const double& alphaSMZin, const int& alphaSorderin, const double& alphaSnfmaxin, const int& iordin);
  void wate96_();
}


// Initialize at the start of the computation
int ReactionRT_DISNC::initAtStart(const string &s)
{
  int isout = Super::initAtStart(s);
  // Basic init:
  
  return isout;
}

// Main function to compute results at an iteration
int ReactionRT_DISNC::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
  return Super::compute(dataSetID,val,err);
}


// 
void ReactionRT_DISNC::initAtIteration() {

  Super::initAtIteration ();
  // optimal:
  double varin[4] = {0.0, 1.0, -2./3., 1.0};

  const double mc = GetParam("mch");
  const double mb = GetParam("mbt");
  const double mZ = GetParam("Mz");
  const double qs0 = 1.0;
  const double as_q0 = alphaS(sqrt(qs0));
  const double as_MZ = alphaS(mZ);
  
  const string order = GetParamS("Order");

  const int  iord = OrderMap( order) - 1;
  const int  asOrederIn = 0;  // ???
  const int  alphaSnfmaxin = 3;

  rt_set_input_(varin, mc, mb, as_q0, as_MZ,  asOrederIn, alphaSnfmaxin, iord);
  wate96_();
 }

// RT
void ReactionRT_DISNC::F2 BASE_PARS 
{
  valarray<double> f2base, f2gamma_base;
  valarray<double> f2gamma_RT(GetNpoint(dataSetID));

  // Get ZMVFNs F2s:
  Super::F2gamma(dataSetID, f2gamma_base, err);
  Super::F2(dataSetID, f2base, err);

  // Get RT F2gamma
  F2gamma_RT(dataSetID, f2gamma_RT, err);

  // Re-scale F2:
  val = f2base * f2gamma_RT / f2gamma_base;

  std::cout << f2gamma_RT[0] << " "<<f2gamma_base[0] << std::endl;
}

void ReactionRT_DISNC::F2gamma_RT BASE_PARS 
{
  // Get x,Q2 arrays:
  auto *q2p  = GetBinValues(dataSetID,"Q2"), *xp  = GetBinValues(dataSetID,"x");
  auto q2 = *q2p, x = *xp;
  
  const size_t Np = GetNpoint(dataSetID);
  int iflag = 1;

  double f2(0), f2b(0), f2c(0), fl(0), flc(0), flb(0);

  for (size_t i=0; i<Np; i++) {
    mstwnc_wrap_(x[i], q2[i], 1,
		 f2, f2c, f2b, fl, flc, flb,
		 iflag, i+1, 1., 0.1, 0 );

    switch ( _dataType) 
      {
      case dataType::sigred :
	val[i] = f2; break;
      case dataType::f2c :
	val[i] = f2c; break ;
      case dataType::f2b :
	val[i] = f2b; break ;
      }
  }
}

