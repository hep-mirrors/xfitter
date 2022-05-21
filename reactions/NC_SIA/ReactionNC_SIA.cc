/*
   @file ReactionFONLL_DISCC.cc
   @date 2017-11-29
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-11-29
*/

#include <iostream>
#include "ReactionNC_SIA.h"
#include "xfitter_pars.h"
// APFEL C++ interface header
#include "APFEL/APFEL.h"
#include "hf_errlog.h"
#include "BaseEvolution.h"
using namespace std;

// The class factory
extern "C" ReactionNC_SIA *create()
{
  return new ReactionNC_SIA();
}

// Initialize at the start of the computation
void ReactionNC_SIA::atStart()
{
  /// ReactionBaseDISCC::atStart();  # this checks QCDNUM
}

void ReactionNC_SIA::initTerm(TermData *td)
{
  ReactionBaseDISCC::initTerm(td);
  _dsIDs[td->id] = td;
  // Check if APFEL evolution is used
  xfitter::BaseEvolution* pdf = td->getPDF();
  if (pdf->getClassName() != string("APFEL") ) {
    std::cerr<<"[ERROR] Reaction "<<getReactionName()<<" only supports APFEL evolution; got evolution named \""<<pdf->_name<<"\" of class \""<<pdf->getClassName()<<"\" for termID="<<td->id<<std::endl;
    hf_errlog(19051815,"F: Reaction "+getReactionName()+" can only work with APFEL evolution, see stderr");
  } 
}

// Compute all predictions in here and store them to be returned
// by the specific functions.
void ReactionNC_SIA::atIteration()
{
  ReactionBaseDISCC::atIteration();
  APFEL::SetProcessDIS("NC");
  // Starting scale ...
  double Q0 = *XFITTER_PARS::getParamD("Q0");

  // Loop over the data sets.
  for (auto tdpair : _dsIDs)
  {
    auto termID = tdpair.first;
    auto td = tdpair.second;
    auto rd = (BaseDISCC::ReactionData *)td->reactionData;
    // Charge of the projectile.
    const double charge = rd->_charge;
    if (charge < 0)
      APFEL::SetProjectileDIS("electron");
    else
      APFEL::SetProjectileDIS("positron");

    // Get x,Q2 arrays.
    APFEL::SetFKObservable("SIA_XSEC");
    auto *q2p = BaseDISCC::GetBinValues(td, "Q2");
    auto *xp = BaseDISCC::GetBinValues(td, "x");
    auto q2 = *q2p;
    auto x = *xp;
    const size_t Np = x.size();
    // Resize arrays.
    _obs[termID].resize(Np);
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
      cout<<"########hello hamed Reaction #############";
        APFEL::ComputeStructureFunctionsAPFEL(Q0, Q);
        APFEL::FKObservables(x[i],Q,0.01);
      }

      // Compute structure functions by interpolation in x for the
      // appropriate component (total, charm, or bottom).

        const double Q = sqrt(q2[i]);
        _obs[termID][i] = APFEL::FKObservables(x[i],Q,0.01);

      Q2save = q2[i];
    }
  }
}

valarray<double> ReactionNC_SIA::OBS(TermData *td)
{
  return _obs[td->id];
}
