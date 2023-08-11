/*
   @file ReactionACOT.cc
   @date 2017-04-10
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-04-10
*/

#include <iostream>
#include "ReactionACOT.h"
#include "xfitter_cpp_base.h"
#include "ext_pdfs.h"

// the class factories
extern "C" ReactionACOT *create()
{
  return new ReactionACOT();
}

// ACOT wrappers from ACOT/src/acotncX_wrap.f:
extern "C"
{
  void acotnc_wrapa_(const double &x, const double &q2, const int &ipn,
		     double &f2, double &f2c, double &f2b, double &fl, double &flc, double &flb);

  void acot_set_input_(const double *varinacot, const int *intvarin );
  
}

// Initialize at the start of the computation
void ReactionACOT::atStart()
{
  Super::atStart();
}

void ReactionACOT::initTerm(TermData *td)
{
  Super::initTerm(td);

  // Allocate internal arrays:
  unsigned termID = td->id;
  _f2acot[termID].resize(GetNpoint(termID));
  _flacot[termID].resize(GetNpoint(termID));
}

//
void ReactionACOT::atIteration()
{
  Super::atIteration();
  // Flag for internal arrays
  for (auto ds : _dsIDs)
  {
    (_f2acot[ds])[0] = -100.;
    (_flacot[ds])[0] = -100.;
  }
}

//
void ReactionACOT::compute(TermData *td, valarray<double> &val, map<string, valarray<double>> &err)
{
  // First init, then call base class:
  td->actualizeWrappers();
  vector<double> varinacot = {*td->getParamD("varinacot0"), *td->getParamD("varinacot1"), *td->getParamD("varinacot2"), *td->getParamD("varinacot3")}; // {0.0, 1.0, -2./3., 1.0};
  vector<int> intvarin = {td->getParamI("intvarin0"), td->getParamI("intvarin1"), td->getParamI("intvarin2"), td->getParamI("intvarin3")}; // {NORD, dum, dum, dum};
  //
  const double mc = *td->getParamD("mch");
  const double mb = *td->getParamD("mbt");
  const double mZ = *td->getParamD("Mz");
  const double qs0 = 1.0;
  const double as_q0 = alphas_wrapper_(sqrt(qs0));
  const double as_MZ = alphas_wrapper_(mZ);

  const string order = td->getParamS("Order");

  const int iord = OrderMap(order) - 1;
  const int asOrederIn = 0; // ???
  const int alphaSnfmaxin = 3;

  // set PDFs, alphaS functions:
  acot_set_pdfs_alphaS(pdf_xfxq_wrapper_, alphas_wrapper_);

  //  acot_set_input_ 
  acot_set_input_(&varinacot[0], &intvarin[0]);


  Super::compute(td, val, err);
}

//
void ReactionACOT::F2 BASE_PARS
{
  unsigned termID = td->id;

  valarray<double> f2base, f2gamma_base;
  valarray<double> f2gamma_ACOT(GetNpoint(termID));

  // Get ACOT F2gamma
  F2gamma_ACOT(td, f2gamma_ACOT, err);

  if (GetDataFlav(termID) == dataFlav::incl)
  {
    // Get ZMVFNs F2s:
    Super::F2gamma(td, f2gamma_base, err);
    Super::F2(td, f2base, err);

    // Re-scale F2:
    val = f2base * f2gamma_ACOT / f2gamma_base;
  }
  else
    val = f2gamma_ACOT;
}

void ReactionACOT::FL BASE_PARS
{
  unsigned termID = td->id;
  valarray<double> flbase, flgamma_base;
  valarray<double> flgamma_ACOT(GetNpoint(termID));

  // Get ACOT F2gamma
  FLgamma_ACOT(td, flgamma_ACOT, err);

  // OZ 19.10.2017 TODO: in dis_sigma.f there is no rescaling for FL at order = 1, should it be here?
  if (GetDataFlav(termID) == dataFlav::incl)
  {
    // Get ZMVFNs F2s:
    Super::FLgamma(td, flgamma_base, err);
    Super::FL(td, flbase, err);

    // Re-scale F2:
    val = flbase * flgamma_ACOT / flgamma_base;
  }
  else
    val = flgamma_ACOT;
}

void ReactionACOT::F2gamma_ACOT BASE_PARS
{
  auto termID = td->id;
  calcF2FL(td);
  val = _f2acot[termID];
}

void ReactionACOT::FLgamma_ACOT BASE_PARS
{
  auto termID = td->id;
  calcF2FL(td);
  val = _flacot[termID];
}

// Place calculations in one function, to optimize calls.
void ReactionACOT::calcF2FL(TermData *td)
{
  unsigned termID = td->id;
  if ((_f2acot[termID][0] < -99.))
  { // compute
    // Get x,Q2 arrays:
    auto *q2p = GetBinValues(td, "Q2"), *xp = GetBinValues(td, "x");
    auto q2 = *q2p, x = *xp;

    const size_t Np = GetNpoint(termID);
    int iflag = 1;

    double f2(0), f2b(0), f2c(0), fl(0), flc(0), flb(0);

    for (size_t i = 0; i < Np; i++)
    {
      if (q2[i] > 1.0)
      {

        acotnc_wrapa_(x[i], q2[i], 1,
		      f2, f2c, f2b, fl, flc, flb);


      }

      switch (GetDataFlav(termID))
      {
      case dataFlav::incl:
        _f2acot[termID][i] = f2;
        _flacot[termID][i] = fl;
        break;
      case dataFlav::c:
        _f2acot[termID][i] = f2c;
        _flacot[termID][i] = flc;
        break;
      case dataFlav::b:
        _f2acot[termID][i] = f2b;
        _flacot[termID][i] = flb;
        break;
      }
    }
  }
}
