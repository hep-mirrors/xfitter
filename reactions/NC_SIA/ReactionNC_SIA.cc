/*
   @file ReactionFONLL_DISCC.cc
   @date 2017-11-29
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-11-29
*/

#include <iostream>
#include "ReactionNC_SIA.h"
#include "xfitter_pars.h"
#include "xfitter_cpp.h"
// APFEL C++ interface header
#include "APFEL/APFEL.h"
#include "hf_errlog.h"
#include "BaseEvolution.h"
#include <cstring>
using namespace std;

// The class factory
extern "C" ReactionNC_SIA *create()
{

  return new ReactionNC_SIA();
}

// Initialize at the start of the computation
void ReactionNC_SIA::atStart()
{
}
void ReactionNC_SIA::initTerm(TermData *td)
{
  unsigned termID = td->id;

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
  //!u/unsigned termID = td->id;
  // dataObs hamed
  _dsIDs.push_back(termID);
  _tdDS[termID] = td;
  _dataObs[termID] = dataObs::incl;

  string msg = "I: calculating SIA cross section";
  //flavObs incl, inclNor, inclNF4, inclNorNF4, inclCH,inclBT
  if (td->hasParam("flavObs"));
  {
    string flavObs = td->getParamS("flavObs");
    if (flavObs == "incl")
    {
      _dataObs[termID] = dataObs::incl;
      msg += "inclusive";
    }
    else if (flavObs == "inclNor")
    {
      _dataObs[termID] = dataObs::inclNor;
      msg += "inclusive_nor";
    }
    else if (flavObs == "inclNF4")
    {
      _dataObs[termID] = dataObs::inclNF4;
      msg += "inclusive_NF4";
    }
    else if (flavObs == "inclNorNF4")
    {
      _dataObs[termID] = dataObs::inclNorNF4;
      msg += "inclusiveNorNf4";
    }
    else if (flavObs == "inclCH")
    {
      _dataObs[termID] = dataObs::inclCH;
      msg += "charm";
    }
    else if (flavObs == "inclBT")
    {
      _dataObs[termID] = dataObs::inclBT;
      msg += "beauty";
    }
    else if (flavObs == "inclLi")
    {
      _dataObs[termID] = dataObs::inclLi;
      msg += "light";
    }
    else
    {
      char buffer[256];
      sprintf(buffer, "F: dataset with id = %d has unknown flavObs = %s", termID, flavObs.c_str());
      string str = buffer;
      hf_errlog(17101902, str);
    }
  }
  auto *q2p = td->getBinColumnOrNull("Q2");
  _npoints[termID] = (*q2p).size();
  _obs[termID].resize(_npoints[termID]);

  // Get PDF id
}
//start hamed
void ReactionNC_SIA::compute(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal)
{

  unsigned termID = td->id;
  valarray<double> val;
  map<string, valarray<double>> err;
  //const double factor = 1.0;
  //val *= factor;
  val = OBS(td);//, val, err);
  valExternal = val;
  errExternal = err;
}
///start
void ReactionNC_SIA::atIteration()
{
  super::atIteration();
  for (auto ds : _dsIDs)
  {
    (_obs[ds])[0] = -100.;
  }

  //ReactionBaseDISNC::atIteration();
  APFEL::SetProcessDIS("NC");
  // Starting scale ...
  double Q0 = *XFITTER_PARS::getParamD("Q0");

  // Loop over the data sets.
  for (auto termID : _dsIDs)
  {
    TermData *td = GetTermData(termID);
   // auto td = tdpair.second;
   // auto rd = (BaseDISCC::ReactionData *)td->reactionData;
    // Non-apfelFF evolution
    if (_non_apfel_evol) {
      td -> actualizeWrappers();
      APFEL::SetPDFSet("external1");
    }
    // Charge of the projectile.
    // This does not have any impact because the difference between
    // "electron" and "positron" for a NC cross section is relevant
    // only when constructing the reduced cross section. But we keep
    // it here for clarity.
    const double charge = GetCharge(termID);    // Charge of the projectile.
//    const double charge = rd->_charge;
    if (charge < 0)
      APFEL::SetProjectileDIS("electron");
    else
      APFEL::SetProjectileDIS("positron");

    // Get x,Q2 arrays.
    // compute the SIA CS by observables
    switch (GetDataObs(termID))
    {
      case dataObs::incl:
          APFEL::SetFKObservable("SIA_XSEC");

          break;
      case dataObs::inclNor:
          APFEL::SetFKObservable("SIA_NORM_XSEC");

          break;
      case dataObs::inclNF4:
          APFEL::SetFKObservable("SIA_XSEC_NF4");

          break;
      case dataObs::inclNorNF4:
          APFEL::SetFKObservable("SIA_NORM_XSEC_NF4");

          break;
      case dataObs::inclCH:
          APFEL::SetFKObservable("SIA_NORM_XSEC_CH");

          break;
      case dataObs::inclBT:
          APFEL::SetFKObservable("SIA_NORM_XSEC_BT");

          break;
      case dataObs::inclLi:
          APFEL::SetFKObservable("SIA_NORM_XSEC_L");

          break;
    }

    //auto *q2p = GetBinValues(td, "Q2");
    //auto *xp = GetBinValues(td, "x");
    auto *q2p = td->getBinColumnOrNull("Q2");
    auto *xp = td->getBinColumnOrNull("x");
    auto q2 = *q2p;
    auto x = *xp;

    const size_t Np = GetNpoint(termID);
    // Resize arrays.
    _obs[termID].resize(Np);
    double Q2save = 0;
    //APFEL::SetFKObservable("SIA_NORM_XSEC");

    for (size_t i = 0; i < Np; i++)
    {
      // Skip all points with Q2 < 1 GeV^2.
      if (q2[i] < 1)
        continue;

      // Recompute structure functions only if the value of Q2
      // changes.
      //if (q2[i] != Q2save)
      //{
        const double Q = sqrt(q2[i]);
        APFEL::ComputeStructureFunctionsAPFEL(Q0, Q);
       // APFEL::FKObservables(x[i],Q,0.01);
      //}

      // Compute structure functions by interpolation in x for the
      // appropriate component (total, charm, or bottom).
       // APFEL::SetFKObservable("SIA_NORM_XSEC");

        //const double Q = sqrt(q2[i]);
        _obs[termID][i] = APFEL::FKObservables(x[i],Q,0.01);
         //cout<<"##hello hamed Reaction 3## "<< x[i]<<" ## "<< Q0 << " ## "<<Q<<"##"<<APFEL::FKObservables(x[i],Q,0.01)<<"***\n";
      Q2save = q2[i];
    }
  }
}
//end hamed


// Compute all predictions in here and store them to be returned
// by the specific functions.



valarray<double> ReactionNC_SIA::OBS(TermData *td)
{
  return _obs[td->id];
  cout<<"##hello hamed Reaction 4##***\n";
}
