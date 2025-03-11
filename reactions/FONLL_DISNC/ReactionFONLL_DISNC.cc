/*
   @file ReactionFONLL_DISNC.cc
   @date 2017-11-29
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-11-29
*/

#include <iostream>
#include "ReactionFONLL_DISNC.h"
#include "xfitter_cpp.h"
#include "xfitter_pars.h"

// APFEL C++ interface header
#include "APFEL/APFEL.h"
#include "hf_errlog.h"
#include "BaseEvolution.h"
#include <cstring>

// The class factory
extern "C" ReactionFONLL_DISNC *create()
{
  return new ReactionFONLL_DISNC();
}

void ReactionFONLL_DISNC::initTerm(TermData* td)
{
  ReactionBaseDISNC::initTerm(td);
  // Check if APFEL evolution is used
  xfitter::BaseEvolution* pdf = td->getPDF();
  if (pdf->getClassName() != string("APFEL") ) {
    hf_errlog(19051815,"I: Reaction " + getReactionName() + " uses non-APFEL evolution. Make sure that APFEL evolution is included");
    _non_apfel_evol = true;

    // We want to check if apfelff was initialized:
    bool foundApfel = false;
    for ( auto const& entry : XFITTER_PARS::gEvolutions ) {
      if (strcmp(entry.second->getClassName(),"APFEL")==0) {
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
void ReactionFONLL_DISNC::atIteration()
{
  ReactionBaseDISNC::atIteration();

  APFEL::SetProcessDIS("NC");

  // Starting scale ...
  double Q0 = *XFITTER_PARS::getParamD("Q0");

  // Loop over the data sets.
  for (auto termID : _dsIDs)
  {
    TermData *td = GetTermData(termID);
    
    // Non-apfelFF evolution
    if (_non_apfel_evol) {
      td -> actualizeWrappers();
      APFEL::SetPDFSet("external1");
    }
    
    // Get evolution boundary:
    double gridLowQLimit = td->getPDF()->getQgrid()[0];
    
    // Charge of the projectile.
    // This does not have any impact because the difference between
    // "electron" and "positron" for a NC cross section is relevant
    // only when constructing the reduced cross section. But we keep
    // it here for clarity.
    const double charge = GetCharge(termID);

    if (charge < 0)
      APFEL::SetProjectileDIS("electron");
    else
      APFEL::SetProjectileDIS("positron");

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
        if (_non_apfel_evol) {
          APFEL::ComputeStructureFunctionsAPFEL(Q, Q); /// do not evolve PDFs, take them from external set at Q
        }
        else {
          APFEL::ComputeStructureFunctionsAPFEL(Q0, Q); /// THIS INCLUDES APFEL EVOLUTION
        }
      }

      // Compute structure functions by interpolation in x for the
      // appropriate component (total, charm, or bottom).
      switch (GetDataFlav(termID))
      {
      case dataFlav::incl:
        _f2fonll[termID][i] = APFEL::F2total(x[i]);
        _flfonll[termID][i] = APFEL::FLtotal(x[i]);
        _f3fonll[termID][i] = -charge * APFEL::F3total(x[i]);
        break;
      case dataFlav::c:
        _f2fonll[termID][i] = APFEL::F2charm(x[i]);
        _flfonll[termID][i] = APFEL::FLcharm(x[i]);
        _f3fonll[termID][i] = -charge * APFEL::F3charm(x[i]);
        break;
      case dataFlav::b:
        _f2fonll[termID][i] = APFEL::F2bottom(x[i]);
        _flfonll[termID][i] = APFEL::FLbottom(x[i]);
        _f3fonll[termID][i] = -charge * APFEL::F3bottom(x[i]);
        break;
      }
      Q2save = q2[i];
    }
  }
}

// FONLL structure functions
void ReactionFONLL_DISNC::F2 BASE_PARS
{
  val = _f2fonll[td->id];
}

void ReactionFONLL_DISNC::FL BASE_PARS
{
  val = _flfonll[td->id];
}

void ReactionFONLL_DISNC::xF3 BASE_PARS
{
  val = _f3fonll[td->id];
}
