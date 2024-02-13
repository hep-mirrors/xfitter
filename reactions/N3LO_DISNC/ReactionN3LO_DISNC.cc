/*
   @file ReactionN3LO_DISNC.cc
*/

#include <iostream>
#include <iomanip>
#include "ReactionN3LO_DISNC.h"
#include "xfitter_pars.h"
#include "xfitter_steer.h"
#include "xfitter_cpp_base.h"
// APFEL C++ interface header
#include <apfel/apfelxx.h>
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
  // x-space grid (the grid parameters should be in parameters.yaml
  const YAML::Node yamlNode=XFITTER_PARS::getEvolutionNode("proton-APFELxx");
  const YAML::Node xGrid = yamlNode["xGrid"];

  vector<apfel::SubGrid> sgv;
  for(auto const& sg : xGrid)
    sgv.push_back(apfel::SubGrid{sg[0].as<int>(), sg[1].as<double>(), sg[2].as<int>()});

  Grid = std::unique_ptr<const apfel::Grid>(new apfel::Grid(sgv));

  // Vectors of thresholds
  const double* MCharm   = XFITTER_PARS::getParamD("mch");
  const double* MBottom  = XFITTER_PARS::getParamD("mbt");
  const double* MTop     = XFITTER_PARS::getParamD("mtp");
  Thresholds = {0, 0, 0, *MCharm, *MBottom, *MTop};

  // Initialize coefficient functions
  F2Obj = InitializeF2NCObjectsZM(*Grid, Thresholds);
  FLObj = InitializeFLNCObjectsZM(*Grid, Thresholds);
  F3Obj = InitializeF3NCObjectsZM(*Grid, Thresholds);
}

void ReactionN3LO_DISNC::initTerm(TermData *td)
{
  ReactionBaseDISNC::initTerm(td);
  // Check if APFELxx evolution is used
  xfitter::BaseEvolution* pdf = td->getPDF();
  if (pdf->getClassName() != string("APFELxx") ) {
    hf_errlog(19051815,"I: Reaction "+getReactionName()+ " uses non-APFELxx evolution. Make sure that APFELxx evolution is included");

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
}

// Compute all predictions in here and store them to be returned
// by the specific functions.
void ReactionN3LO_DISNC::atIteration()
{
  ReactionBaseDISNC::atIteration();

  // Perturbative order
  const int PerturbativeOrder    = OrderMap(XFITTER_PARS::getParamS("Order")) - 1;

  // Effective charges
  std::function<std::vector<double>(double const&)> fBq = [=] (double const& Q) -> std::vector<double> { return apfel::ElectroWeakCharges(Q, false); };
  std::function<std::vector<double>(double const&)> fDq = [=] (double const& Q) -> std::vector<double> { return apfel::ParityViolatingElectroWeakCharges(Q, false); };

  // Evolved PDFs and alphas from BaseEvolution
  xfitter::BaseEvolution* pdf = xfitter::get_evolution();
  const auto PDFs = [&] (double const& x, double const& Q) -> std::map<int, double> { return apfel::PhysToQCDEv(pdf->xfxQmap(x, Q)); };
  const auto as = [&] (double const& Q) -> double { return pdf->getAlphaS(Q); };

  // Initialize structure functions
  const auto F2 = BuildStructureFunctions(F2Obj, PDFs, PerturbativeOrder, as, fBq);
  const auto FL = BuildStructureFunctions(FLObj, PDFs, PerturbativeOrder, as, fBq);
  const auto F3 = BuildStructureFunctions(F3Obj, PDFs, PerturbativeOrder, as, fDq);

  //ZM+M-M0
    
  // Tabulate Structure functions
  const apfel::TabulateObject<apfel::Distribution> F2total {[&] (double const& Q) -> apfel::Distribution{ return F2.at(0).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> FLtotal {[&] (double const& Q) -> apfel::Distribution{ return FL.at(0).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F3total {[&] (double const& Q) -> apfel::Distribution{ return F3.at(0).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};

  // Loop over the data sets.
  for (auto termID : _dsIDs)
    {
      TermData *td = GetTermData(termID);

      if (GetDataFlav(termID) != dataFlav::incl)
	continue;
    
      // Get evolution boundary:
      double gridLowQLimit = td->getPDF()->getQgrid()[0];

      // Charge of the projectile.
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

      for (size_t i = 0; i < Np; i++)
	{
	  // Skip all points with Q2 below starting scale
	  if (q2[i] < gridLowQLimit*gridLowQLimit)
	    continue;
	
	  // Compute structure functions by interpolation in x and Q
	  _f2fonll[termID][i] = F2total.EvaluatexQ(x[i], sqrt(q2[i]));
	  _flfonll[termID][i] = FLtotal.EvaluatexQ(x[i], sqrt(q2[i]));
	  _f3fonll[termID][i] = -charge * F3total.EvaluatexQ(x[i], sqrt(q2[i]));
	}
    }

  bool initcharm = false;
  for (auto termID : _dsIDs)
    if (GetDataFlav(termID) == dataFlav::c)
      initcharm = true;

  if (initcharm)
    {
      const apfel::TabulateObject<apfel::Distribution> F2charm {[&] (double const& Q) -> apfel::Distribution{ return F2.at(4).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> FLcharm {[&] (double const& Q) -> apfel::Distribution{ return FL.at(4).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> F3charm {[&] (double const& Q) -> apfel::Distribution{ return F3.at(4).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};

      // Loop over the data sets.
      for (auto termID : _dsIDs)
	{
	  TermData *td = GetTermData(termID);

	  if (GetDataFlav(termID) != dataFlav::c)
	    continue;
    
	  // Get evolution boundary:
	  double gridLowQLimit = td->getPDF()->getQgrid()[0];

	  // Charge of the projectile.
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

	  for (size_t i = 0; i < Np; i++)
	    {
	      // Skip all points with Q2 below starting scale
	      if (q2[i] < gridLowQLimit*gridLowQLimit)
		continue;
	
	      // Compute structure functions by interpolation in x and Q
	      _f2fonll[termID][i] = F2charm.EvaluatexQ(x[i], sqrt(q2[i]));
	      _flfonll[termID][i] = FLcharm.EvaluatexQ(x[i], sqrt(q2[i]));
	      _f3fonll[termID][i] = -charge * F3charm.EvaluatexQ(x[i], sqrt(q2[i]));
	    }
	}
    }
  
  
  bool initbottom = false;
  for (auto termID : _dsIDs)
    if (GetDataFlav(termID) == dataFlav::b)
      initbottom = true;

  if (initbottom)
    {
      const apfel::TabulateObject<apfel::Distribution> F2bottom{[&] (double const& Q) -> apfel::Distribution{ return F2.at(5).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> FLbottom{[&] (double const& Q) -> apfel::Distribution{ return FL.at(5).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> F3bottom{[&] (double const& Q) -> apfel::Distribution{ return F3.at(5).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};

      // Loop over the data sets.
      for (auto termID : _dsIDs)
	{
	  TermData *td = GetTermData(termID);

	  if (GetDataFlav(termID) != dataFlav::b)
	    continue;
    
	  // Get evolution boundary:
	  double gridLowQLimit = td->getPDF()->getQgrid()[0];

	  // Charge of the projectile.
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

	  for (size_t i = 0; i < Np; i++)
	    {
	      // Skip all points with Q2 below starting scale
	      if (q2[i] < gridLowQLimit*gridLowQLimit)
		continue;
	
	      // Compute structure functions by interpolation in x and Q
	      _f2fonll[termID][i] = F2bottom.EvaluatexQ(x[i], sqrt(q2[i]));
	      _flfonll[termID][i] = FLbottom.EvaluatexQ(x[i], sqrt(q2[i]));
	      _f3fonll[termID][i] = -charge * F3bottom.EvaluatexQ(x[i], sqrt(q2[i]));
	    }
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
