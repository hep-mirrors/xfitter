/*
   @file ReactionN3LO_DISNC.cc
*/

#include <iostream>
#include <iomanip>
#include "ReactionN3LO_DISNC.h"
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
extern "C" ReactionN3LO_DISNC *create()
{
  return new ReactionN3LO_DISNC();
}

// Initialize at the start of the computation
void ReactionN3LO_DISNC::atStart()
{
  /// ReactionBaseDISNC::atStart();  # this checks QCDNUM

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100,1e-5,3}, apfel::SubGrid{60,1e-1,3}, apfel::SubGrid{50,6e-1,3}, apfel::SubGrid{50,8e-1,3}}};

  // Vectors of thresholds
  const double* MCharm   = XFITTER_PARS::getParamD("mch");
  const double* MBottom  = XFITTER_PARS::getParamD("mbt");
  const double* MTop     = XFITTER_PARS::getParamD("mtp");
  //const std::vector<double> Thresholds = {0, 0, 0, *MCharm, *MBottom, *MTop};
  Thresholds = {0, 0, 0, *MCharm, *MBottom, *MTop};

  // Initialize coefficient functions
  F2Obj = InitializeF2NCObjectsZM(g, Thresholds);
  FLObj = InitializeFLNCObjectsZM(g, Thresholds);
  F3Obj = InitializeF3NCObjectsZM(g, Thresholds);
}

void ReactionN3LO_DISNC::initTerm(TermData *td)
{
  ReactionBaseDISNC::initTerm(td);
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
void ReactionN3LO_DISNC::atIteration()
{
  ReactionBaseDISNC::atIteration();

  // Starting scale ...
  double Q0 = *XFITTER_PARS::getParamD("Q0");

  // Perturbative order
  const int PerturbativeOrder    = OrderMap(XFITTER_PARS::getParamS("Order")) - 1;

  // Effective charges
  std::function<std::vector<double>(double const&)> fBq = [=] (double const& Q) -> std::vector<double> { return apfel::ElectroWeakCharges(Q, false); };
  std::function<std::vector<double>(double const&)> fDq = [=] (double const& Q) -> std::vector<double> { return apfel::ParityViolatingElectroWeakCharges(Q, false); };

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100,1e-5,3}, apfel::SubGrid{60,1e-1,3}, apfel::SubGrid{50,6e-1,3}, apfel::SubGrid{50,8e-1,3}}};
  
  // Loop over the data sets.
  for (auto termID : _dsIDs)
  {
    TermData *td = GetTermData(termID);

    // Non-apfelFF evolution
    if (_non_apfel_evol) {
      td -> actualizeWrappers();
      //APFEL::SetPDFSet("external1");
    }

    // Initialize coefficient functions
    const auto F2Obj = InitializeF2NCObjectsZM(g, Thresholds);
    const auto FLObj = InitializeFLNCObjectsZM(g, Thresholds);
    const auto F3Obj = InitializeF3NCObjectsZM(g, Thresholds);
    
    xfitter::EvolutionAPFELxx* pdf = (xfitter::EvolutionAPFELxx*) td->getPDF();
    const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs = pdf->GetTabulatedPDFs();
    const std::function<double(double const& Q)> as = pdf->GetAlphaQCD();

    // Evolved PDFs
    const auto PDFs = [&] (double const& x, double const& Q) -> std::map<int, double> { return TabulatedPDFs.EvaluateMapxQ(x, Q); };

    // Initialize structure functions
    const auto F2 = BuildStructureFunctions(F2Obj, PDFs, PerturbativeOrder, as, fBq);
    const auto FL = BuildStructureFunctions(FLObj, PDFs, PerturbativeOrder, as, fBq);
    const auto F3 = BuildStructureFunctions(F3Obj, PDFs, PerturbativeOrder, as, fDq);

    // Tabulate Structure functions
    const apfel::TabulateObject<apfel::Distribution> F2total {[&] (double const& Q) -> apfel::Distribution{ return F2.at(0).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
    const apfel::TabulateObject<apfel::Distribution> F2light {[&] (double const& Q) -> apfel::Distribution{ return F2.at(1).Evaluate(Q) + F2.at(2).Evaluate(Q) + F2.at(3).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
    const apfel::TabulateObject<apfel::Distribution> F2charm {[&] (double const& Q) -> apfel::Distribution{ return F2.at(4).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
    const apfel::TabulateObject<apfel::Distribution> F2bottom{[&] (double const& Q) -> apfel::Distribution{ return F2.at(5).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};

    const apfel::TabulateObject<apfel::Distribution> FLtotal {[&] (double const& Q) -> apfel::Distribution{ return FL.at(0).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
    const apfel::TabulateObject<apfel::Distribution> FLlight {[&] (double const& Q) -> apfel::Distribution{ return FL.at(1).Evaluate(Q) + FL.at(2).Evaluate(Q) + FL.at(3).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
    const apfel::TabulateObject<apfel::Distribution> FLcharm {[&] (double const& Q) -> apfel::Distribution{ return FL.at(4).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
    const apfel::TabulateObject<apfel::Distribution> FLbottom{[&] (double const& Q) -> apfel::Distribution{ return FL.at(5).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};

    const apfel::TabulateObject<apfel::Distribution> F3total {[&] (double const& Q) -> apfel::Distribution{ return F3.at(0).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
    const apfel::TabulateObject<apfel::Distribution> F3light {[&] (double const& Q) -> apfel::Distribution{ return F3.at(1).Evaluate(Q) + F3.at(2).Evaluate(Q) + F3.at(3).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
    const apfel::TabulateObject<apfel::Distribution> F3charm {[&] (double const& Q) -> apfel::Distribution{ return F3.at(4).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
    const apfel::TabulateObject<apfel::Distribution> F3bottom{[&] (double const& Q) -> apfel::Distribution{ return F3.at(5).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};

    // Get evolution boundary:
    double gridLowQLimit = td->getPDF()->getQgrid()[0];

    // Charge of the projectile.
    // This does not have any impact because the difference between
    // "electron" and "positron" for a NC cross section is relevant
    // only when constructing the reduced cross section. But we keep
    // it here for clarity.
    const double charge = GetCharge(termID);

    // Get x,Q2 arrays.
    auto *q2p = GetBinValues(td, "Q2");
    auto *xp = GetBinValues(td, "x");
    auto q2 = *q2p;
    auto x = *xp;

    const size_t Np = GetNpoint(termID);
    // Resize arrays.
    _f2fonll[termID].resize(Np);
    _flfonll[termID].resize(Np);
    _f3fonll[termID].resize(Np);

    double Q2save = 0;
    for (size_t i = 0; i < Np; i++)
    {
      // Skip all points with Q2 below starting scale
      if (q2[i] < gridLowQLimit*gridLowQLimit)
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
      switch (GetDataFlav(termID))
      {
      case dataFlav::incl:
        _f2fonll[termID][i] = F2total.EvaluatexQ(x[i], sqrt(q2[i]));// / 2;
        _flfonll[termID][i] = FLtotal.EvaluatexQ(x[i], sqrt(q2[i]));// / 2;
        _f3fonll[termID][i] = -charge * F3total.EvaluatexQ(x[i], sqrt(q2[i]));// / 2;
	break;
      case dataFlav::c:
        _f2fonll[termID][i] = F2charm.EvaluatexQ(x[i], sqrt(q2[i]));// / 2;
        _flfonll[termID][i] = FLcharm.EvaluatexQ(x[i], sqrt(q2[i]));// / 2;
        _f3fonll[termID][i] = -charge * F3charm.EvaluatexQ(x[i], sqrt(q2[i]));// / 2;
        break;
      case dataFlav::b:
        _f2fonll[termID][i] = F2bottom.EvaluatexQ(x[i], sqrt(q2[i]));
        _flfonll[termID][i] = FLbottom.EvaluatexQ(x[i], sqrt(q2[i]));
        _f3fonll[termID][i] = -charge * F3bottom.EvaluatexQ(x[i], sqrt(q2[i]));
        break;
      }

      Q2save = q2[i];
    }
  }
}

// N3LO structure functions
void ReactionN3LO_DISNC::F2 BASE_PARS
{
  val = _f2fonll[td->id];
}

void ReactionN3LO_DISNC::FL BASE_PARS
{
  val = _flfonll[td->id];
}

void ReactionN3LO_DISNC::xF3 BASE_PARS
{
  val = _f3fonll[td->id];
}
