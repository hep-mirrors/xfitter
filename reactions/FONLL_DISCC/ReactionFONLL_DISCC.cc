/*
   @file ReactionFONLL_DISCC.cc
   @date 2017-11-29
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-11-29
*/

#include <iostream>
#include "ReactionFONLL_DISCC.h"
#include "xfitter_pars.h"
// APFEL C++ interface header
#include "APFEL/APFEL.h"
#include "hf_errlog.h"
#include "BaseEvolution.h"

// The class factory
extern "C" ReactionFONLL_DISCC *create()
{
  return new ReactionFONLL_DISCC();
}

// Initialize at the start of the computation
void ReactionFONLL_DISCC::atStart()
{
  /// ReactionBaseDISCC::atStart();  # this checks QCDNUM
}

void ReactionFONLL_DISCC::initTerm(TermData *td)
{
  ReactionBaseDISCC::initTerm(td);
  // Also prepare a map of DS:
  _dsIDs[td->id] = td;
  // Check if APFEL evolution is used
  xfitter::BaseEvolution* pdf = td->getPDF();
  if (pdf->getClassName() != string("APFEL") ) {
    hf_errlog(19051815,"I: Reaction "+getReactionName()+ " uses non-APFEL evolution. Make sure that APFEL evolution is included");
    _non_apfel_evol = true;

    // We want to check if apfelff was initialized:
    bool foundApfel = false;
    for ( auto const& entry : XFITTER_PARS::gEvolutions ) {
      if (entry.second->getClassName() == string("APFEL")) {
	foundApfel = true;
	break;
      }
    }
    if (not foundApfel) {
      hf_errlog(21122901,"F: Include fortran APFEL evolution to Evolutions: list (even if it is not used) to initialize FONLL DIS modules");
    }
  }
  else {
    _non_apfel_evol = false;
  }    
}

// Compute all predictions in here and store them to be returned
// by the specific functions.
void ReactionFONLL_DISCC::atIteration()
{
  ReactionBaseDISCC::atIteration();
  APFEL::SetProcessDIS("CC");

  // Starting scale ...
  double Q0 = *XFITTER_PARS::getParamD("Q0");

  // Loop over the data sets.
  for (auto tdpair : _dsIDs)
  {
    auto termID = tdpair.first;
    auto td = tdpair.second;

    // Non-apfelFF evolution
    if (_non_apfel_evol) {
      td -> actualizeWrappers();
      APFEL::SetPDFSet("external1");
    }
        
    auto rd = (BaseDISCC::ReactionData *)td->reactionData;
    // Charge of the projectile.
    const double charge = rd->_charge;
    if (charge < 0)
      APFEL::SetProjectileDIS("electron");
    else
      APFEL::SetProjectileDIS("positron");

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
        APFEL::ComputeStructureFunctionsAPFEL(Q, Q); /// do not evolve PDFs, take them from external set at Q
      }

      // Compute structure functions by interpolation in x for the
      // appropriate component (total, charm, or bottom).
      switch (rd->_dataFlav)
      {
      case BaseDISCC::dataFlav::incl:
        _f2fonll[termID][i] = APFEL::F2total(x[i]) / 2;
        _flfonll[termID][i] = APFEL::FLtotal(x[i]) / 2;
        _f3fonll[termID][i] = APFEL::F3total(x[i]) / 2;
        break;
      case BaseDISCC::dataFlav::c:
        _f2fonll[termID][i] = APFEL::F2charm(x[i]) / 2;
        _flfonll[termID][i] = APFEL::FLcharm(x[i]) / 2;
        _f3fonll[termID][i] = APFEL::F3charm(x[i]) / 2;
        break;
      }

      Q2save = q2[i];
    }
  }
}

// FONLL structure functions
valarray<double> ReactionFONLL_DISCC::F2(TermData *td)
{
  return _f2fonll[td->id];
}

valarray<double> ReactionFONLL_DISCC::FL(TermData *td)
{
  return _flfonll[td->id];
}

valarray<double> ReactionFONLL_DISCC::xF3(TermData *td)
{
  return _f3fonll[td->id];
}
