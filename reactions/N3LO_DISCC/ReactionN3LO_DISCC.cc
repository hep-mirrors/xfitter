/*
   @file ReactionN3LO_DISCC.cc
*/

#include <iostream>
#include <iomanip>
#include "ReactionN3LO_DISCC.h"
#include <cmath>
#include <algorithm>
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
  apfel::SetVerbosityLevel(0);

  // x-space grid
  const YAML::Node Node     = XFITTER_PARS::rootNode["byReaction"];
  const YAML::Node yamlNode = Node["N3LO_DISCC"];
  const YAML::Node xGrid = yamlNode["xGrid"];

  vector<apfel::SubGrid> sgv;
  for(auto const& sg : xGrid)
    sgv.push_back(apfel::SubGrid{sg[0].as<int>(), sg[1].as<double>(), sg[2].as<int>()});

  Grid = std::unique_ptr<const apfel::Grid>(new apfel::Grid(sgv));
  _grid = Grid.get();

  // FK mode, as N3LO_DISNC
  if (yamlNode["FK"]) _fk = yamlNode["FK"].as<int>() != 0;
  if (yamlNode["FKcheck"]) _fkCheck = yamlNode["FKcheck"].as<int>();
  if (_fk) {
    xfitter::BaseEvolution* evolution = xfitter::get_evolution();
    if (evolution->getClassName() != string("APFELxx"))
      hf_errlog(26102306, "F: Reaction N3LO_DISCC: FK mode needs the APFELxx evolution as default evolution");
    _apfelxx = static_cast<xfitter::EvolutionAPFELxx*>(evolution);
    const YAML::Node evolutionGrid = XFITTER_PARS::getEvolutionNode(evolution->_name)["xGrid"];
    bool same = evolutionGrid.size() == xGrid.size();
    for (size_t i = 0; same && i < xGrid.size(); ++i)
      same = evolutionGrid[i][0].as<int>() == xGrid[i][0].as<int>() && evolutionGrid[i][1].as<double>() == xGrid[i][1].as<double>()
             && evolutionGrid[i][2].as<int>() == xGrid[i][2].as<int>();
    if (!same)
      hf_errlog(26102307, "F: Reaction N3LO_DISCC: FK mode needs the same xGrid as the APFELxx evolution");
    _grid = &_apfelxx->GetGrid();
    _apfelxx->RequestOperators();
    hf_errlog(26102308, "I: Reaction N3LO_DISCC: FK mode, structure functions from FK tables of the APFELxx evolution operators");
  }

  // Vectors of thresholds
  const double* MCharm   = XFITTER_PARS::getParamD("mch");
  const double* MBottom  = XFITTER_PARS::getParamD("mbt");
  const double* MTop     = XFITTER_PARS::getParamD("mtp");
  Thresholds = {0, 0, 0, *MCharm, *MBottom, *MTop};

  // Initialize coefficient functions
  F2PlusCCObj  = InitializeF2CCPlusObjectsZM (*_grid, Thresholds);
  F2MinusCCObj = InitializeF2CCMinusObjectsZM(*_grid, Thresholds);
  FLPlusCCObj  = InitializeFLCCPlusObjectsZM (*_grid, Thresholds);
  FLMinusCCObj = InitializeFLCCMinusObjectsZM(*_grid, Thresholds);
  F3PlusCCObj  = InitializeF3CCPlusObjectsZM (*_grid, Thresholds);
  F3MinusCCObj = InitializeF3CCMinusObjectsZM(*_grid, Thresholds);
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
  if (!_fk) {
    computeFromPDFs();
    return;
  }
  if (_fkVersion != _apfelxx->OperatorVersion()) {
    buildFK();
    _fkVersion = _apfelxx->OperatorVersion();
    if (_fkCheck) {
      computeFromPDFs();
      const auto f2 = _f2fonll, fl = _flfonll, f3 = _f3fonll;
      computeFromFK();
      double worst = 0;
      for (auto const& t : _f2fonll)
        for (size_t i = 0; i < t.second.size(); ++i)
          for (auto const& pair : {std::make_pair(f2.at(t.first)[i], _f2fonll[t.first][i]),
                                   std::make_pair(fl.at(t.first)[i], _flfonll[t.first][i]),
                                   std::make_pair(f3.at(t.first)[i], _f3fonll[t.first][i])})
            if (pair.first != 0) worst = std::max(worst, std::abs(pair.second / pair.first - 1));
      std::cout << "N3LO_DISCC FK check: max |F_FK / F - 1| = " << worst << std::endl;
      if (!(worst < 1e-4))
        hf_errlog(26102309, "W: Reaction N3LO_DISCC: FK tables differ from the default computation by more than 1e-4");
      return;
    }
  }
  computeFromFK();
}

// Default computation: structure functions from the evolved PDFs
void ReactionN3LO_DISCC::computeFromPDFs()
{
  bool tabulate = true;

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

	  //non tabulated version
	  if (!tabulate)
	    {
	      _f2fonll[termID][i] =          F2p.at(0).Evaluate(x[i], sqrt(q2[i])) + charge*F2m.at(0).Evaluate(x[i], sqrt(q2[i]));
	      _flfonll[termID][i] =          FLp.at(0).Evaluate(x[i], sqrt(q2[i])) + charge*FLm.at(0).Evaluate(x[i], sqrt(q2[i]));
	      _f3fonll[termID][i] = charge * F3p.at(0).Evaluate(x[i], sqrt(q2[i])) +        F3m.at(0).Evaluate(x[i], sqrt(q2[i]));
	    }
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
// FK mode: rows of F2, FL, xF3 per point (FKTables.h), as computeFromPDFs
void ReactionN3LO_DISCC::buildFK()
{
  using xfitter::fk::OperatorTable;
  using xfitter::fk::TabulateFK;
  const YAML::Node QGrid = XFITTER_PARS::rootNode["byReaction"]["N3LO_DISCC"]["QGrid"];
  const int n = QGrid[0].as<int>();
  const double qmin = QGrid[1].as<double>(), qmax = QGrid[2].as<double>();
  const int ord = QGrid[3].as<int>();
  const int PO = OrderMap(XFITTER_PARS::getParamS("Order")) - 1;
  const std::function<std::vector<double>(double const&)> fCKM = [] (double const&) -> std::vector<double> { return apfel::CKM2; };
  xfitter::EvolutionAPFELxx* e = _apfelxx;
  const std::function<double(double const&)> as = [e] (double const& Q) { return e->getAlphaS(Q); };
  OperatorTable const& E5 = _apfelxx->GetOperatorTable(5);
  const auto tab = [&] (decltype(F2PlusCCObj) const& obj, std::vector<int> const& k) {
    return TabulateFK(obj, k, PO, as, fCKM, E5, n, qmin, qmax, ord, Thresholds);
  };
  const int nJoint = static_cast<int>(_apfelxx->InputDistributions().at(0).GetDistributionJointGrid().size());
  _fkRows.clear();
  for (const bool charmOnly : {false, true}) {
    const auto flav = charmOnly ? BaseDISCC::dataFlav::c : BaseDISCC::dataFlav::incl;
    bool any = false;
    for (auto tdpair : _dsIDs)
      any = any || ((BaseDISCC::ReactionData*)tdpair.second->reactionData)->_dataFlav == flav;
    if (!any) continue;
    const std::vector<int> k = charmOnly ? std::vector<int>{4, 5} : std::vector<int>{0};
    const OperatorTable F2p = tab(F2PlusCCObj, k), F2m = tab(F2MinusCCObj, k);
    const OperatorTable FLp = tab(FLPlusCCObj, k), FLm = tab(FLMinusCCObj, k);
    const OperatorTable F3p = tab(F3PlusCCObj, k), F3m = tab(F3MinusCCObj, k);
    for (auto tdpair : _dsIDs) {
      auto td = tdpair.second;
      auto rd = (BaseDISCC::ReactionData*)td->reactionData;
      if (rd->_dataFlav != flav) continue;
      const double charge = rd->_charge;
      const auto q2 = *BaseDISCC::GetBinValues(td, "Q2");
      const auto x = *BaseDISCC::GetBinValues(td, "x");
      const size_t Np = x.size();
      FKRows& r = _fkRows[tdpair.first];
      r.f2.assign(Np, xfitter::fk::Row(nJoint));
      r.fl.assign(Np, xfitter::fk::Row(nJoint));
      r.f3.assign(Np, xfitter::fk::Row(nJoint));
      r.active.assign(Np, false);
      for (size_t i = 0; i < Np; ++i) {
        if (q2[i] < 1) continue;
        r.active[i] = true;
        const double Q = sqrt(q2[i]);
        r.f2[i].add(F2p, x[i], Q, 1);       r.f2[i].add(F2m, x[i], Q, charge);
        r.fl[i].add(FLp, x[i], Q, 1);       r.fl[i].add(FLm, x[i], Q, charge);
        r.f3[i].add(F3p, x[i], Q, charge);  r.f3[i].add(F3m, x[i], Q, 1);
      }
    }
  }
}

void ReactionN3LO_DISCC::computeFromFK()
{
  const std::vector<double> f0 = xfitter::fk::Flatten(_apfelxx->InputDistributions());
  for (auto const& t : _fkRows) {
    const FKRows& r = t.second;
    const size_t Np = r.active.size();
    _f2fonll[t.first].resize(Np);
    _flfonll[t.first].resize(Np);
    _f3fonll[t.first].resize(Np);
    for (size_t i = 0; i < Np; ++i) {
      _f2fonll[t.first][i] = r.active[i] ? r.f2[i].contract(f0) : 0.;
      _flfonll[t.first][i] = r.active[i] ? r.fl[i].contract(f0) : 0.;
      _f3fonll[t.first][i] = r.active[i] ? r.f3[i].contract(f0) : 0.;
    }
  }
}

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
