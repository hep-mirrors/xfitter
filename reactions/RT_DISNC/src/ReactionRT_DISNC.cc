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
  // use external instead of include:
  void rt_set_pdfs_alphaS( pXFXlike xfx, pOneParFunc aS);    //!< Set PDFs and alphaS
}


// Initialize at the start of the computation
int ReactionRT_DISNC::atStart(const string &s)
{
  int isout = Super::atStart(s);
  // Basic init:

  return isout;
}

void ReactionRT_DISNC::setDatasetParameters( int dataSetID, map<string,string> pars, map<string,double> parsDataset) {
  Super::setDatasetParameters(dataSetID, pars, parsDataset);
  // Allocate internal arrays:
  _f2rt[dataSetID].resize(GetNpoint(dataSetID));
  _flrt[dataSetID].resize(GetNpoint(dataSetID));
}

//
void ReactionRT_DISNC::initAtIteration() {

  Super::initAtIteration ();
  // optimal:
  vector<double> varin = GetParamV("varin"); // {0.0, 1.0, -2./3., 1.0};
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

  // set PDFs, alphaS functions:
  rt_set_pdfs_alphaS( getXFX(), getAlphaS() );

  rt_set_input_(&varin[0], mc, mb, as_q0, as_MZ,  asOrederIn, alphaSnfmaxin, iord);
  wate96_();
  // Flag for internal arrays
  for ( auto ds : _dsIDs)  {
    (_f2rt[ds])[0] = -100.;
    (_flrt[ds])[0] = -100.;
  }

 }

// RT
void ReactionRT_DISNC::F2 BASE_PARS
{
  valarray<double> f2base, f2gamma_base;
  valarray<double> f2gamma_RT(GetNpoint(dataSetID));

  // Get RT F2gamma
  F2gamma_RT(dataSetID, f2gamma_RT, err);

  if(GetDataFlav(dataSetID) == dataFlav::incl)
  {
    // Get ZMVFNs F2s:
    Super::F2gamma(dataSetID, f2gamma_base, err);
    Super::F2(dataSetID, f2base, err);

    // Re-scale F2:
    val = f2base * f2gamma_RT / f2gamma_base;
  }
  else
    val = f2gamma_RT;
}

void ReactionRT_DISNC::FL BASE_PARS
{
  valarray<double> flbase, flgamma_base;
  valarray<double> flgamma_RT(GetNpoint(dataSetID));

  // Get RT F2gamma
  FLgamma_RT(dataSetID, flgamma_RT, err);

  // OZ 19.10.2017 TODO: in dis_sigma.f there is no rescaling for FL at order = 1, should it be here?
  if(GetDataFlav(dataSetID) == dataFlav::incl)
  {
    // Get ZMVFNs F2s:
    Super::FLgamma(dataSetID, flgamma_base, err);
    Super::FL(dataSetID, flbase, err);

    // Re-scale F2:
    val = flbase * flgamma_RT / flgamma_base;
  }
  else
    val = flgamma_RT;
}

void ReactionRT_DISNC::F2gamma_RT BASE_PARS
{
  calcF2FL(dataSetID);
  val = _f2rt[dataSetID];
}

void ReactionRT_DISNC::FLgamma_RT BASE_PARS
{
  calcF2FL(dataSetID);
  val = _flrt[dataSetID];
}


// Place calculations in one function, to optimize calls.
void ReactionRT_DISNC::calcF2FL(int dataSetID) {
  if ( (_f2rt[dataSetID][0]< -99.) ) { // compute
  // Get x,Q2 arrays:
    auto *q2p  = GetBinValues(dataSetID,"Q2"), *xp  = GetBinValues(dataSetID,"x");
    auto q2 = *q2p, x = *xp;

    const size_t Np = GetNpoint(dataSetID);
    int iflag = 1;

    double f2(0), f2b(0), f2c(0), fl(0), flc(0), flb(0);

    for (size_t i=0; i<Np; i++) {
        if (q2[i]>1.0) {

            mstwnc_wrap_(x[i], q2[i], 1,
                         f2, f2c, f2b, fl, flc, flb,
                         iflag, i+1, 1., 0.1, 0 );
        }


      switch ( GetDataFlav(dataSetID) )
      {
      case dataFlav::incl :
	_f2rt[dataSetID][i] = f2;
	_flrt[dataSetID][i] = fl;
	break;
      case dataFlav::c :
	_f2rt[dataSetID][i] = f2c;
	_flrt[dataSetID][i] = flc;
	break ;
      case dataFlav::b :
	_f2rt[dataSetID][i] = f2b;
	_flrt[dataSetID][i] = flb;
	break ;
      }
    }
  }
}
