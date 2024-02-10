/*
   @file ReactionN3LO_DISCC.cc
*/

#include <iostream>
#include <iomanip>
#include "ReactionN3LO_DISCC.h"
#include "xfitter_pars.h"
#include "xfitter_cpp_base.h"
// APFEL C++ interface header
#include <apfel/apfelxx.h>
//#include <apfel/alphaqcd.h>
//#include <apfel/messages.h>
//#include <apfel/rotations.h>
#include "hf_errlog.h"
#include "BaseEvolution.h"
#include "EvolutionAPFELxx.h"
#include "EvolutionQCDNUM.h"

// The class factory
extern "C" ReactionN3LO_DISCC *create()
{
  return new ReactionN3LO_DISCC();
}

// Initialize at the start of the computation
void ReactionN3LO_DISCC::atStart()
{
  /// ReactionBaseDISCC::atStart();  # this checks QCDNUM

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100,1e-5,3}, apfel::SubGrid{60,1e-1,3}, apfel::SubGrid{50,6e-1,3}, apfel::SubGrid{50,8e-1,3}}};

  // Vectors of thresholds
  const double* MCharm   = XFITTER_PARS::getParamD("mch");
  const double* MBottom  = XFITTER_PARS::getParamD("mbt");
  const double* MTop     = XFITTER_PARS::getParamD("mtp");
  //const std::vector<double> Thresholds = {0, 0, 0, *MCharm, *MBottom, *MTop};
  Thresholds = {0, 0, 0, *MCharm, *MBottom, *MTop};

  // Initialize coefficient functions
  F2PlusCCObj  = InitializeF2CCPlusObjectsZM(g, Thresholds);
  F2MinusCCObj = InitializeF2CCMinusObjectsZM(g, Thresholds);
  FLPlusCCObj  = InitializeFLCCPlusObjectsZM(g, Thresholds);
  FLMinusCCObj = InitializeFLCCMinusObjectsZM(g, Thresholds);
  F3PlusCCObj  = InitializeF3CCPlusObjectsZM(g, Thresholds);
  F3MinusCCObj = InitializeF3CCMinusObjectsZM(g, Thresholds);
}

void ReactionN3LO_DISCC::initTerm(TermData *td)
{
  ReactionBaseDISCC::initTerm(td);
  // Also prepare a map of DS:
  _dsIDs[td->id] = td;
  // Check if APFELxx evolution is used
  xfitter::BaseEvolution* pdf = td->getPDF();
  if (pdf->getClassName() != string("APFELxx") ) {
    hf_errlog(19051815,"I: Reaction "+getReactionName()+ " uses non-APFELxx evolution. Make sure that APFELxx evolution is included");
    _non_apfel_evol = true;

    // We want to check if apfelxx was initialized:
    bool foundApfel = false;
    for ( auto const& entry : XFITTER_PARS::gEvolutions ) {
      if (entry.second->getClassName() == string("APFELxx")) {
	foundApfel = true;
	break;
      }
    }
    if (not foundApfel) {
      hf_errlog(21122901,"F: Include APFELxx evolution to Evolutions: list (even if it is not used) to initialize N3LO DIS modules");
    }
  }
  else {
    _non_apfel_evol = false;
  }    
}

// Compute all predictions in here and store them to be returned
// by the specific functions.
void ReactionN3LO_DISCC::atIteration()
{
  ReactionBaseDISCC::atIteration();

  // Starting scale ...
  double Q0 = *XFITTER_PARS::getParamD("Q0");


  // CKM matrix elements
  std::function<std::vector<double>(double const&)> fCKM = [=] (double const&) -> std::vector<double> { return apfel::CKM2; };

  // Initial scale
  const double mu0 = Q0;

  // Perturbative order
  const int PerturbativeOrder    = OrderMap(XFITTER_PARS::getParamS("Order")) - 1;

//  // Running coupling
//  apfel::AlphaQCD a{0.35, sqrt(2), Thresholds, PerturbativeOrder};
//  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
//  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  
  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100,1e-5,3}, apfel::SubGrid{60,1e-1,3}, apfel::SubGrid{50,6e-1,3}, apfel::SubGrid{50,8e-1,3}}};
  
  // Loop over the data sets.
  for (auto tdpair : _dsIDs)
  {
    auto termID = tdpair.first;
    auto td = tdpair.second;

    // Non-apfelFF evolution
    if (_non_apfel_evol) {
      td -> actualizeWrappers();
      //APFEL::SetPDFSet("external1");
    }

    // Initialize coefficient functions
    const auto F2PlusCCObj  = InitializeF2CCPlusObjectsZM(g, Thresholds);
    const auto F2MinusCCObj = InitializeF2CCMinusObjectsZM(g, Thresholds);
    const auto FLPlusCCObj  = InitializeFLCCPlusObjectsZM(g, Thresholds);
    const auto FLMinusCCObj = InitializeFLCCMinusObjectsZM(g, Thresholds);
    const auto F3PlusCCObj  = InitializeF3CCPlusObjectsZM(g, Thresholds);
    const auto F3MinusCCObj = InitializeF3CCMinusObjectsZM(g, Thresholds);
    
    xfitter::EvolutionAPFELxx* pdf = (xfitter::EvolutionAPFELxx*) td->getPDF();
    const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs = pdf->GetTabulatedPDFs();
    const std::function<double(double const& Q)> as = pdf->GetAlphaQCD();

    // // Initialize QCD evolution objects
    // const auto DglapObj = InitializeDglapObjectsQCD(g, Thresholds);
    // 
    // // Construct the DGLAP object
    // auto EvolvedPDFs = BuildDglap(DglapObj, apfel::LHToyPDFs, mu0, PerturbativeOrder, as);
    // 
    // // Tabulate PDFs
    // const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3};
    
    // Evolved PDFs
    const auto PDFs = [&] (double const& x, double const& Q) -> std::map<int, double> { return TabulatedPDFs.EvaluateMapxQ(x, Q); };

    // Initialize structure functions
    const auto F2p = BuildStructureFunctions(F2PlusCCObj,  PDFs, PerturbativeOrder, as, fCKM);
    const auto F2m = BuildStructureFunctions(F2MinusCCObj, PDFs, PerturbativeOrder, as, fCKM);
    const auto FLp = BuildStructureFunctions(FLPlusCCObj,  PDFs, PerturbativeOrder, as, fCKM);
    const auto FLm = BuildStructureFunctions(FLMinusCCObj, PDFs, PerturbativeOrder, as, fCKM);
    const auto F3p = BuildStructureFunctions(F3PlusCCObj,  PDFs, PerturbativeOrder, as, fCKM);
    const auto F3m = BuildStructureFunctions(F3MinusCCObj, PDFs, PerturbativeOrder, as, fCKM);

    const apfel::TabulateObject<apfel::Distribution> F2total  {[&] (double const& Q) -> apfel::Distribution { return F2p.at(0).Evaluate(Q) - F2m.at(0).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
    const apfel::TabulateObject<apfel::Distribution> F2charm  {[&] (double const& Q) -> apfel::Distribution { return F2p.at(4).Evaluate(Q) - F2m.at(4).Evaluate(Q) + F2p.at(5).Evaluate(Q) - F2m.at(5).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
    const apfel::TabulateObject<apfel::Distribution> FLtotal  {[&] (double const& Q) -> apfel::Distribution { return FLp.at(0).Evaluate(Q) - FLm.at(0).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
    const apfel::TabulateObject<apfel::Distribution> FLcharm  {[&] (double const& Q) -> apfel::Distribution { return FLp.at(4).Evaluate(Q) - FLm.at(4).Evaluate(Q) + FLp.at(5).Evaluate(Q) - FLm.at(5).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
    const apfel::TabulateObject<apfel::Distribution> F3total  {[&] (double const& Q) -> apfel::Distribution { return F3p.at(0).Evaluate(Q) - F3m.at(0).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
    const apfel::TabulateObject<apfel::Distribution> F3charm  {[&] (double const& Q) -> apfel::Distribution { return F3p.at(4).Evaluate(Q) - F3m.at(4).Evaluate(Q) + F3p.at(5).Evaluate(Q) - F3m.at(5).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};

    //const apfel::TabulateObject<apfel::Distribution> F2light  {[&] (double const& Q) -> apfel::Distribution { return F2p.at(1).Evaluate(Q) - F2m.at(1).Evaluate(Q) + F2p.at(2).Evaluate(Q) - F2m.at(2).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
    //const apfel::TabulateObject<apfel::Distribution> F2bottom {[&] (double const& Q) -> apfel::Distribution { return F2p.at(3).Evaluate(Q) - F2m.at(3).Evaluate(Q) + F2p.at(6).Evaluate(Q) - F2m.at(6).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
    //const apfel::TabulateObject<apfel::Distribution> FLlight  {[&] (double const& Q) -> apfel::Distribution { return FLp.at(1).Evaluate(Q) - FLm.at(1).Evaluate(Q) + FLp.at(2).Evaluate(Q) - FLm.at(2).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
    //const apfel::TabulateObject<apfel::Distribution> FLbottom {[&] (double const& Q) -> apfel::Distribution { return FLp.at(3).Evaluate(Q) - FLm.at(3).Evaluate(Q) + FLp.at(6).Evaluate(Q) - FLm.at(6).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
    //const apfel::TabulateObject<apfel::Distribution> F3light  {[&] (double const& Q) -> apfel::Distribution { return F3p.at(1).Evaluate(Q) - F3m.at(1).Evaluate(Q) + F3p.at(2).Evaluate(Q) - F3m.at(2).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
    //const apfel::TabulateObject<apfel::Distribution> F3bottom {[&] (double const& Q) -> apfel::Distribution { return F3p.at(3).Evaluate(Q) - F3m.at(3).Evaluate(Q) + F3p.at(6).Evaluate(Q) - F3m.at(6).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
    
    auto rd = (BaseDISCC::ReactionData *)td->reactionData;
    // Charge of the projectile.
    const double charge = rd->_charge;
//    if (charge < 0)
//      APFEL::SetProjectileDIS("electron");
//    else
//      APFEL::SetProjectileDIS("positron");

    // Get x,Q2 arrays.
    auto *q2p = BaseDISCC::GetBinValues(td, "Q2");
    auto *xp = BaseDISCC::GetBinValues(td, "x");
    auto q2 = *q2p;
    auto x = *xp;

    const size_t Np = x.size();
    // Resize arrays.
    _f2fonll[termID].resize(Np);
    _flfonll[termID].resize(Np);
    _f3fonll[termID].resize(Np);

    double Q2save = 0;
    for (size_t i = 0; i < Np; i++)
    {
      // Skip all points with Q2 < 1 GeV^2.
      if (q2[i] < 1)
        continue;

      // Recompute structure functions only if the value of Q2
      // changes.
      if (q2[i] != Q2save)
      {
        const double Q = sqrt(q2[i]);
        //APFEL::ComputeStructureFunctionsAPFEL(Q0, Q);
      }

      // Compute structure functions by interpolation in x for the
      // appropriate component (total, charm, or bottom).
      switch (rd->_dataFlav)
      {
      case BaseDISCC::dataFlav::incl:
        _f2fonll[termID][i] = F2total.EvaluatexQ(x[i], sqrt(q2[i]));// / 2;
        _flfonll[termID][i] = FLtotal.EvaluatexQ(x[i], sqrt(q2[i]));// / 2;
        _f3fonll[termID][i] = F3total.EvaluatexQ(x[i], sqrt(q2[i]));// / 2;

	break;
      case BaseDISCC::dataFlav::c:
        _f2fonll[termID][i] = F2charm.EvaluatexQ(x[i], sqrt(q2[i]));// / 2;
        _flfonll[termID][i] = FLcharm.EvaluatexQ(x[i], sqrt(q2[i]));// / 2;
        _f3fonll[termID][i] = F3charm.EvaluatexQ(x[i], sqrt(q2[i]));// / 2;
        break;
      }

      Q2save = q2[i];
    }
  }
}

// N3LO structure functions
valarray<double> ReactionN3LO_DISCC::F2(TermData *td)
{
  return _f2fonll[td->id];
}

valarray<double> ReactionN3LO_DISCC::FL(TermData *td)
{
  return _flfonll[td->id];
}

valarray<double> ReactionN3LO_DISCC::xF3(TermData *td)
{
  return _f3fonll[td->id];
}
