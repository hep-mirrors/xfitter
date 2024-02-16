/*
   @file ReactionN3LO_DISCC.cc
*/

#include <iostream>
#include <iomanip>
#include "ReactionN3LO_DISCC.h"
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
extern "C" ReactionN3LO_DISCC *create()
{
  return new ReactionN3LO_DISCC();
}

// Initialize at the start of the computation
void ReactionN3LO_DISCC::atStart()
{
  // x-space grid
  const YAML::Node Node     = XFITTER_PARS::rootNode["byReaction"];
  const YAML::Node yamlNode = Node["N3LO_DISCC"];
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
  F2PlusCCObj  = InitializeF2CCPlusObjectsZM (*Grid, Thresholds);
  F2MinusCCObj = InitializeF2CCMinusObjectsZM(*Grid, Thresholds);
  FLPlusCCObj  = InitializeFLCCPlusObjectsZM (*Grid, Thresholds);
  FLMinusCCObj = InitializeFLCCMinusObjectsZM(*Grid, Thresholds);
  F3PlusCCObj  = InitializeF3CCPlusObjectsZM (*Grid, Thresholds);
  F3MinusCCObj = InitializeF3CCMinusObjectsZM(*Grid, Thresholds);
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
void ReactionN3LO_DISCC::atIteration()
{
  ReactionBaseDISCC::atIteration();

  // CKM matrix elements --> should read from yaml
  const double Vud = *(XFITTER_PARS::getParamD("Vud"));
  const double Vus = *(XFITTER_PARS::getParamD("Vus"));
  const double Vub = *(XFITTER_PARS::getParamD("Vub"));
  const double Vcd = *(XFITTER_PARS::getParamD("Vcd"));
  const double Vcs = *(XFITTER_PARS::getParamD("Vcs"));
  const double Vcb = *(XFITTER_PARS::getParamD("Vcb"));
  const double Vtd = *(XFITTER_PARS::getParamD("Vtd"));
  const double Vts = *(XFITTER_PARS::getParamD("Vts"));
  const double Vtb = *(XFITTER_PARS::getParamD("Vtb"));
  const std::vector<double> CKM2 = {pow(Vud,2), pow(Vus,2), pow(Vub,2), pow(Vcd,2), pow(Vcs,2), pow(Vcb,2), pow(Vtd,2), pow(Vts,2), pow(Vtb,2)};
  std::function<std::vector<double>(double const&)> fCKM = [=] (double const&) -> std::vector<double> { return apfel::CKM2; };

  //Q grid parameters
  const YAML::Node Node     = XFITTER_PARS::rootNode["byReaction"];
  const YAML::Node yamlNode = Node["N3LO_DISCC"];
  const YAML::Node QGrid    = yamlNode["QGrid"];
  
  int n = QGrid[0].as<int>();
  double qmin = QGrid[1].as<double>();
  double qmax = QGrid[2].as<double>();
  int ord = QGrid[3].as<int>();
  
  // Perturbative order
  const int PerturbativeOrder    = OrderMap(XFITTER_PARS::getParamS("Order")) - 1;

  // Evolved PDFs and alphas from BaseEvolution
  xfitter::BaseEvolution* pdf = xfitter::get_evolution();
  const auto PDFs = [&] (double const& x, double const& Q) -> std::map<int, double> { return apfel::PhysToQCDEv(pdf->xfxQmap(x, Q)); };
  const auto as = [&] (double const& Q) -> double { return pdf->getAlphaS(Q); };

  // Initialize structure functions
  const auto F2p = BuildStructureFunctions(F2PlusCCObj,  PDFs, PerturbativeOrder, as, fCKM);
  const auto F2m = BuildStructureFunctions(F2MinusCCObj, PDFs, PerturbativeOrder, as, fCKM);
  const auto FLp = BuildStructureFunctions(FLPlusCCObj,  PDFs, PerturbativeOrder, as, fCKM);
  const auto FLm = BuildStructureFunctions(FLMinusCCObj, PDFs, PerturbativeOrder, as, fCKM);
  const auto F3p = BuildStructureFunctions(F3PlusCCObj,  PDFs, PerturbativeOrder, as, fCKM);
  const auto F3m = BuildStructureFunctions(F3MinusCCObj, PDFs, PerturbativeOrder, as, fCKM);

  const apfel::TabulateObject<apfel::Distribution> F2totalp  {[&] (double const& Q) -> apfel::Distribution { return F2p.at(0).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> FLtotalp  {[&] (double const& Q) -> apfel::Distribution { return FLp.at(0).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F3totalp  {[&] (double const& Q) -> apfel::Distribution { return F3p.at(0).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F2totalm  {[&] (double const& Q) -> apfel::Distribution { return F2m.at(0).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> FLtotalm  {[&] (double const& Q) -> apfel::Distribution { return FLm.at(0).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F3totalm  {[&] (double const& Q) -> apfel::Distribution { return F3m.at(0).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};

  // Loop over the data sets.
  for (auto tdpair : _dsIDs)
    {
      auto termID = tdpair.first;
      auto td = tdpair.second;
      auto rd = (BaseDISCC::ReactionData *)td->reactionData;
      if (rd->_dataFlav != BaseDISCC::dataFlav::incl)
	continue;
    
      // Charge of the projectile.
      const double charge = rd->_charge;

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

      for (size_t i = 0; i < Np; i++)
	{
	  // Skip all points with Q2 < 1 GeV^2.
	  if (q2[i] < 1)
	    continue;

	  // Compute structure functions by interpolation in x and Q
	  _f2fonll[termID][i] =          F2totalp.EvaluatexQ(x[i], sqrt(q2[i])) + charge*F2totalm.EvaluatexQ(x[i], sqrt(q2[i]));
	  _flfonll[termID][i] =          FLtotalp.EvaluatexQ(x[i], sqrt(q2[i])) + charge*FLtotalm.EvaluatexQ(x[i], sqrt(q2[i]));
	  _f3fonll[termID][i] = charge * F3totalp.EvaluatexQ(x[i], sqrt(q2[i])) +        F3totalm.EvaluatexQ(x[i], sqrt(q2[i]));
	}
    }

  bool initcharm = false;
  for (auto tdpair : _dsIDs)
  {
    auto td = tdpair.second;
    auto rd = (BaseDISCC::ReactionData *)td->reactionData;
    if (rd->_dataFlav == BaseDISCC::dataFlav::c)
      initcharm = true;
  }

  if (initcharm)
    {
      const apfel::TabulateObject<apfel::Distribution> F2charmp  {[&] (double const& Q) -> apfel::Distribution { return F2p.at(4).Evaluate(Q) + F2p.at(5).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> FLcharmp  {[&] (double const& Q) -> apfel::Distribution { return FLp.at(4).Evaluate(Q) + FLp.at(5).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> F3charmp  {[&] (double const& Q) -> apfel::Distribution { return F3p.at(4).Evaluate(Q) + F3p.at(5).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> F2charmm  {[&] (double const& Q) -> apfel::Distribution { return F2m.at(4).Evaluate(Q) + F2m.at(5).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> FLcharmm  {[&] (double const& Q) -> apfel::Distribution { return FLm.at(4).Evaluate(Q) + FLm.at(5).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> F3charmm  {[&] (double const& Q) -> apfel::Distribution { return F3m.at(4).Evaluate(Q) + F3m.at(5).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};

      // Loop over the data sets.
      for (auto tdpair : _dsIDs)
	{
	  auto termID = tdpair.first;
	  auto td = tdpair.second;
	  auto rd = (BaseDISCC::ReactionData *)td->reactionData;
	  if (rd->_dataFlav != BaseDISCC::dataFlav::c)
	    continue;
    
	  // Charge of the projectile.
	  const double charge = rd->_charge;

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

	  for (size_t i = 0; i < Np; i++)
	    {
	      // Skip all points with Q2 < 1 GeV^2.
	      if (q2[i] < 1)
		continue;

	      // Compute structure functions by interpolation in x and Q
	      _f2fonll[termID][i] =          F2charmp.EvaluatexQ(x[i], sqrt(q2[i])) + charge*F2charmm.EvaluatexQ(x[i], sqrt(q2[i]));
	      _flfonll[termID][i] =          FLcharmp.EvaluatexQ(x[i], sqrt(q2[i])) + charge*FLcharmm.EvaluatexQ(x[i], sqrt(q2[i]));
	      _f3fonll[termID][i] = charge * F3charmp.EvaluatexQ(x[i], sqrt(q2[i])) +        F3charmm.EvaluatexQ(x[i], sqrt(q2[i]));
	    }
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
