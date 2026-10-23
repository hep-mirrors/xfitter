/*
   @file ReactionN3LO_DISNC.cc
*/

#include <iostream>
#include <iomanip>
#include "ReactionN3LO_DISNC.h"
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

using namespace std;

// The class factory
extern "C" ReactionN3LO_DISNC *create()
{
  return new ReactionN3LO_DISNC();
}

// Initialize at the start of the computation
void ReactionN3LO_DISNC::atStart()
{
  apfel::SetVerbosityLevel(0);
  
  // x-space grid
  const YAML::Node Node     = XFITTER_PARS::rootNode["byReaction"];
  const YAML::Node yamlNode = Node["N3LO_DISNC"];
  const YAML::Node xGrid = yamlNode["xGrid"];

  vector<apfel::SubGrid> sgv;
  for(auto const& sg : xGrid)
    sgv.push_back(apfel::SubGrid{sg[0].as<int>(), sg[1].as<double>(), sg[2].as<int>()});

  Grid = std::unique_ptr<const apfel::Grid>(new apfel::Grid(sgv));
  _grid = Grid.get();

  // FK mode: the rows are built from the APFELxx evolution operators, which need the
  // structure-function objects on the evolution grid
  if (yamlNode["FK"]) _fk = yamlNode["FK"].as<int>() != 0;
  if (yamlNode["FKcheck"]) _fkCheck = yamlNode["FKcheck"].as<int>();
  if (_fk) {
    xfitter::BaseEvolution* evolution = xfitter::get_evolution();
    if (evolution->getClassName() != string("APFELxx"))
      hf_errlog(26102302, "F: Reaction N3LO_DISNC: FK mode needs the APFELxx evolution as default evolution");
    _apfelxx = static_cast<xfitter::EvolutionAPFELxx*>(evolution);
    const YAML::Node evolutionGrid = XFITTER_PARS::getEvolutionNode(evolution->_name)["xGrid"];
    bool same = evolutionGrid.size() == xGrid.size();
    for (size_t i = 0; same && i < xGrid.size(); ++i)
      same = evolutionGrid[i][0].as<int>() == xGrid[i][0].as<int>() && evolutionGrid[i][1].as<double>() == xGrid[i][1].as<double>()
             && evolutionGrid[i][2].as<int>() == xGrid[i][2].as<int>();
    if (!same)
      hf_errlog(26102303, "F: Reaction N3LO_DISNC: FK mode needs the same xGrid as the APFELxx evolution");
    _grid = &_apfelxx->GetGrid();
    _apfelxx->RequestOperators();
    hf_errlog(26102304, "I: Reaction N3LO_DISNC: FK mode, structure functions from FK tables of the APFELxx evolution operators");
  }

  // Vectors of thresholds (5F, no top) and masses
  const double* MCharm   = XFITTER_PARS::getParamD("mch");
  const double* MBottom  = XFITTER_PARS::getParamD("mbt");
  const double* MTop     = XFITTER_PARS::getParamD("mtp");
  Thresholds = {0, 0, 0, *MCharm, *MBottom};
  Masses     = Thresholds;

  // Thresholds of the 3F and 4F PDFs and alphas (for tabulation)
  Thresholds3F = {0, 0, 0};
  Thresholds4F = {0, 0, 0, *MCharm};

  // Masses for the massive coefficient functions: APFELxx requires all
  // 6 masses, the number of zero masses sets the number of light flavours
  Masses3F = {0, 0, 0, *MCharm, *MBottom, *MTop};
  Masses4F = {0, 0, 0, 0,       *MBottom, *MTop};

  // Initialize coefficient functions
  F2Obj = InitializeF2NCObjectsZM(*_grid, Thresholds);
  FLObj = InitializeFLNCObjectsZM(*_grid, Thresholds);
  F3Obj = InitializeF3NCObjectsZM(*_grid, Thresholds);

  F2thr = yamlNode["F2thr"].as<int>();
  F2asy = yamlNode["F2asy"].as<int>();
  F2dmp = yamlNode["F2dmp"].as<int>();
  FLthr = yamlNode["FLthr"].as<int>();
  FLasy = yamlNode["FLasy"].as<int>();
  FLdmp = yamlNode["FLdmp"].as<int>();

  // 3-flavour scheme massive coefficient functions
  F2ObjM3F  = InitializeF2NCObjectsMassive(*_grid,     Masses3F, 1e-5, 150, 0.05, 10000, 3, 0.0005, F2thr, F2asy, F2dmp);
  FLObjM3F  = InitializeFLNCObjectsMassive(*_grid,     Masses3F, 1e-5, 150, 0.05, 10000, 3, 0.0005, FLthr, FLasy, FLdmp);
  F2ObjM03F = InitializeF2NCObjectsMassiveZero(*_grid, Masses3F);
  FLObjM03F = InitializeFLNCObjectsMassiveZero(*_grid, Masses3F);

  // 4-flavour scheme massive coefficient functions
  F2ObjM4F  = InitializeF2NCObjectsMassive(*_grid,     Masses4F, 1e-5, 150, 0.05, 10000, 3, 0.0005, F2thr, F2asy, F2dmp);
  FLObjM4F  = InitializeFLNCObjectsMassive(*_grid,     Masses4F, 1e-5, 150, 0.05, 10000, 3, 0.0005, FLthr, FLasy, FLdmp);
  F2ObjM04F = InitializeF2NCObjectsMassiveZero(*_grid, Masses4F);
  FLObjM04F = InitializeFLNCObjectsMassiveZero(*_grid, Masses4F);
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
  if (!_fk) {
    computeFromPDFs();
    return;
  }
  bool check = false;
  if (_fkVersion != _apfelxx->OperatorVersion()) {
    buildFK();
    _fkVersion = _apfelxx->OperatorVersion();
    check = _fkCheck != 0;
  }
  if (check) {
    // the default computation once, to check the FK tables
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
    std::cout << "N3LO_DISNC FK check: max |F_FK / F - 1| = " << worst << std::endl;
    if (!(worst < 1e-4))
      hf_errlog(26102305, "W: Reaction N3LO_DISNC: FK tables differ from the default computation by more than 1e-4");
    return;
  }
  computeFromFK();
}

// Default computation: structure functions from the evolved PDFs
void ReactionN3LO_DISNC::computeFromPDFs()
{
  bool tabulate = true;
  bool massive = true;

  //Q grid parameters
  const YAML::Node Node     = XFITTER_PARS::rootNode["byReaction"];
  const YAML::Node yamlNode = Node["N3LO_DISNC"];
  const YAML::Node QGrid    = yamlNode["QGrid"];

	// SZ 2024.09.06 for zero-mass scheme: read "massive" from yaml
	if (yamlNode["massive"]) 
	{
		massive = yamlNode["massive"].as<int>();
	}

  int n = QGrid[0].as<int>();
  double qmin = QGrid[1].as<double>();
  double qmax = QGrid[2].as<double>();
  int ord = QGrid[3].as<int>();

  //put in atConfiguration()
  //if (XFITTER_PARS::gParameters.count("N3LO_DISNC/F2thr") != F2thr ||
  //    XFITTER_PARS::gParameters.count("N3LO_DISNC/F2asy") != F2asy ||
  //    XFITTER_PARS::gParameters.count("N3LO_DISNC/F2dmp") != F2dmp)
  //  {
  //    F2thr = *(XFITTER_PARS::gParameters.at("N3LO_DISNC/F2thr"));
  //    F2asy = *(XFITTER_PARS::gParameters.at("N3LO_DISNC/F2asy"));
  //    F2dmp = *(XFITTER_PARS::gParameters.at("N3LO_DISNC/F2dmp"));
  //    F2ObjM3F  = InitializeF2NCObjectsMassive(*_grid,     Masses3F, 1e-5, 150, 0.05, 10000, 3, 0.0005, F2thr, F2asy, F2dmp);
  //  }
  //
  //if (XFITTER_PARS::gParameters.count("N3LO_DISNC/FLthr") != FLthr ||
  //    XFITTER_PARS::gParameters.count("N3LO_DISNC/FLasy") != FLasy ||
  //    XFITTER_PARS::gParameters.count("N3LO_DISNC/FLdmp") != FLdmp)
  //  {
  //    FLthr = *(XFITTER_PARS::gParameters.at("N3LO_DISNC/FLthr"));
  //    FLasy = *(XFITTER_PARS::gParameters.at("N3LO_DISNC/FLasy"));
  //    FLdmp = *(XFITTER_PARS::gParameters.at("N3LO_DISNC/FLdmp"));
  //    FLObjM3F  = InitializeFLNCObjectsMassive(*_grid,     Masses3F, 1e-5, 150, 0.05, 10000, 3, 0.0005, FLthr, FLasy, FLdmp);
  //  }

  // Perturbative order
  const int PerturbativeOrder    = OrderMap(XFITTER_PARS::getParamS("Order")) - 1;

  // Effective charges --> could read from yaml
  //const double sin2thw = *(XFITTER_PARS::getParamD("sin2thW"));
  //const double gf = *(XFITTER_PARS::getParamD("gf"));
  //const double Mz = *(XFITTER_PARS::getParamD("Mz"));
  //const double Mw = *(XFITTER_PARS::getParamD("Mw"));
  //const double Wz = *(XFITTER_PARS::getParamD("Wz"));
  //const double Ww = *(XFITTER_PARS::getParamD("Ww"));
  std::function<std::vector<double>(double const&)> fBq = [=] (double const& Q) -> std::vector<double> { return apfel::ElectroWeakCharges(Q, false); };
  std::function<std::vector<double>(double const&)> fDq = [=] (double const& Q) -> std::vector<double> { return apfel::ParityViolatingElectroWeakCharges(Q, false); };

  // Evolved PDFs and alphas from BaseEvolution
  xfitter::BaseEvolution* pdf = xfitter::get_evolution();
  const auto PDFs = [&] (double const& x, double const& Q) -> std::map<int, double> { return apfel::PhysToQCDEv(pdf->xfxQmap(x, Q)); };
  const auto as = [&] (double const& Q) -> double { return pdf->getAlphaS(Q); };

  // 3F and 4F PDFs and alphas for the massive structure functions
  // (only available with APFELxx evolution, otherwise fall back to the VFNS ones)
  xfitter::EvolutionAPFELxx* apfelxx = nullptr;
  if (pdf->getClassName() == string("APFELxx"))
    apfelxx = static_cast<xfitter::EvolutionAPFELxx*>(pdf);
  else
    hf_errlog(26101303,"W: Reaction "+getReactionName()+": default evolution is not APFELxx, 3F and 4F PDFs and alphas not available, use VFNS PDFs and alphas in massive structure functions");

  const auto PDFs3F = [&] (double const& x, double const& Q) -> std::map<int, double> { return apfel::PhysToQCDEv(apfelxx ? apfelxx->xfxQmapNF(3, x, Q) : pdf->xfxQmap(x, Q)); };
  const auto as3F = [&] (double const& Q) -> double { return apfelxx ? apfelxx->getAlphaSNF(3, Q) : pdf->getAlphaS(Q); };

  const auto PDFs4F = [&] (double const& x, double const& Q) -> std::map<int, double> { return apfel::PhysToQCDEv(apfelxx ? apfelxx->xfxQmapNF(4, x, Q) : pdf->xfxQmap(x, Q)); };
  const auto as4F = [&] (double const& Q) -> double { return apfelxx ? apfelxx->getAlphaSNF(4, Q) : pdf->getAlphaS(Q); };
  
  // Initialize structure functions
  const auto F2 = BuildStructureFunctions(F2Obj, PDFs, PerturbativeOrder, as, fBq);
  const auto FL = BuildStructureFunctions(FLObj, PDFs, PerturbativeOrder, as, fBq);
  const auto F3 = BuildStructureFunctions(F3Obj, PDFs, PerturbativeOrder, as, fDq);
  
  if (!massive)
    {
      // Tabulate Structure functions
      const apfel::TabulateObject<apfel::Distribution> F2total {[&] (double const& Q) -> apfel::Distribution{ return F2.at(0).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> FLtotal {[&] (double const& Q) -> apfel::Distribution{ return FL.at(0).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> F3total {[&] (double const& Q) -> apfel::Distribution{ return F3.at(0).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};

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
    }
  else //general mass VFN -> FONLL
    {
      int msvordmax = 3;
      if (yamlNode["msvordmax"]) 
	{
	  msvordmax = yamlNode["msvordmax"].as<int>();
	}
      // Initialize structure functions
      const auto F2M3F  = BuildStructureFunctions(F2ObjM3F,  PDFs3F, std::min(PerturbativeOrder,msvordmax), as3F, fBq);
      const auto FLM3F  = BuildStructureFunctions(FLObjM3F,  PDFs3F, std::min(PerturbativeOrder,msvordmax), as3F, fBq);
      const auto F2M03F = BuildStructureFunctions(F2ObjM03F, PDFs3F, std::min(PerturbativeOrder,msvordmax), as3F, fBq);
      const auto FLM03F = BuildStructureFunctions(FLObjM03F, PDFs3F, std::min(PerturbativeOrder,msvordmax), as3F, fBq);

      const auto F2M4F  = BuildStructureFunctions(F2ObjM4F,  PDFs4F, std::min(PerturbativeOrder,msvordmax), as4F, fBq);
      const auto FLM4F  = BuildStructureFunctions(FLObjM4F,  PDFs4F, std::min(PerturbativeOrder,msvordmax), as4F, fBq);
      const auto F2M04F = BuildStructureFunctions(F2ObjM04F, PDFs4F, std::min(PerturbativeOrder,msvordmax), as4F, fBq);
      const auto FLM04F = BuildStructureFunctions(FLObjM04F, PDFs4F, std::min(PerturbativeOrder,msvordmax), as4F, fBq);
      
      // Tabulate Structure functions
      const apfel::TabulateObject<apfel::Distribution> F2light {[&] (double const& Q) -> apfel::Distribution{ return F2.at(1).Evaluate(Q) + F2.at(2).Evaluate(Q) + F2.at(3).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> FLlight {[&] (double const& Q) -> apfel::Distribution{ return FL.at(1).Evaluate(Q) + FL.at(2).Evaluate(Q) + FL.at(3).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> F3light {[&] (double const& Q) -> apfel::Distribution{ return F3.at(1).Evaluate(Q) + F3.at(2).Evaluate(Q) + F3.at(3).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};
      
      const apfel::TabulateObject<apfel::Distribution> F2charmZM {[&] (double const& Q) -> apfel::Distribution{ return F2.at(4).Evaluate(Q);   }, n, qmin, qmax, ord, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> F2charmM  {[&] (double const& Q) -> apfel::Distribution{ return F2M3F.at(4).Evaluate(Q);  }, n, qmin, qmax, ord, Thresholds3F};
      const apfel::TabulateObject<apfel::Distribution> F2charmM0 {[&] (double const& Q) -> apfel::Distribution{ return F2M03F.at(4).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds3F};
      const apfel::TabulateObject<apfel::Distribution> FLcharmZM {[&] (double const& Q) -> apfel::Distribution{ return FL.at(4).Evaluate(Q);   }, n, qmin, qmax, ord, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> FLcharmM  {[&] (double const& Q) -> apfel::Distribution{ return FLM3F.at(4).Evaluate(Q);  }, n, qmin, qmax, ord, Thresholds3F};
      const apfel::TabulateObject<apfel::Distribution> FLcharmM0 {[&] (double const& Q) -> apfel::Distribution{ return FLM03F.at(4).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds3F};
      const apfel::TabulateObject<apfel::Distribution> F3charmZM {[&] (double const& Q) -> apfel::Distribution{ return F3.at(4).Evaluate(Q);   }, n, qmin, qmax, ord, Thresholds};
      
      const apfel::TabulateObject<apfel::Distribution> F2bottomZM{[&] (double const& Q) -> apfel::Distribution{ return F2.at(5).Evaluate(Q);   }, n, qmin, qmax, ord, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> F2bottomM {[&] (double const& Q) -> apfel::Distribution{ return F2M4F.at(5).Evaluate(Q);  }, n, qmin, qmax, ord, Thresholds4F};
      const apfel::TabulateObject<apfel::Distribution> F2bottomM0{[&] (double const& Q) -> apfel::Distribution{ return F2M04F.at(5).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds4F};
      const apfel::TabulateObject<apfel::Distribution> FLbottomZM{[&] (double const& Q) -> apfel::Distribution{ return FL.at(5).Evaluate(Q);   }, n, qmin, qmax, ord, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> FLbottomM {[&] (double const& Q) -> apfel::Distribution{ return FLM4F.at(5).Evaluate(Q);  }, n, qmin, qmax, ord, Thresholds4F};
      const apfel::TabulateObject<apfel::Distribution> FLbottomM0{[&] (double const& Q) -> apfel::Distribution{ return FLM04F.at(5).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds4F};
      const apfel::TabulateObject<apfel::Distribution> F3bottomZM{[&] (double const& Q) -> apfel::Distribution{ return F3.at(5).Evaluate(Q);   }, n, qmin, qmax, ord, Thresholds};

      //double xx = 0.0001;
      //string ff = "f2";
      //string qq = "c";
      //int imod = yamlNode["imod"].as<int>();
      //string tag = "";
      //if (imod == 1 || imod == 3 || imod == 5)
      //	tag = "d";
      //if (imod == 2 || imod == 4 || imod == 6)
      //	tag = "u";
      //cout << endl;
      //cout << "{" << endl;
      //cout << "TGraph *g" << ff << qq << PerturbativeOrder << tag << " = new TGraph();" << endl;
      //double Q, logQ, logQ_min = 0, logQ_max = 2;
      //int N = 100;
      //double dlogQ = (logQ_max - logQ_min) / N;
      //for (int i = 0; i < N; i++)
      //	{
      //	  logQ = logQ_min + i * dlogQ;
      //	  Q = pow(10, logQ);
      //	  double etab = Q*Q / ( Q*Q + 4 * Masses[4] * Masses[4] );
      //	  cout << "g" << ff << qq << PerturbativeOrder << tag << "->SetPoint(g" << ff << qq << PerturbativeOrder << tag << "->GetN(), " << Q << ", " << F2charmM.EvaluatexQ(xx/etab, Q) << ");" << endl;
      //	  //	cout << "g" << ff << qq << PerturbativeOrder << tag << "->SetPoint(g" << ff << qq << PerturbativeOrder << tag << "->GetN(), " << Q << ", " << F2M.at(5).Evaluate(xx, Q)   << ");" << endl;
      //	}
      //cout << endl;
      //exit(0);
    
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
	      //ZM+switch(Q/m)*(M-M0)
	      double mc2 = pow(Thresholds[3],2);
	      double fthrc = q2[i] > mc2 ? pow(1.-mc2/q2[i],2): 0.;
	      double mb2 = pow(Thresholds[4],2);
	      double fthrb = q2[i] > mb2 ? pow(1.-mb2/q2[i],2): 0.;
	      double etac = q2[i] / ( q2[i] + 4 * Masses[3] * Masses[3] );
	      double etab = q2[i] / ( q2[i] + 4 * Masses[4] * Masses[4] );
	      
	      _f2fonll[termID][i] =           (+ F2light.EvaluatexQ(x[i], sqrt(q2[i]))
	      				       + fthrc*F2charmZM.EvaluatexQ(x[i], sqrt(q2[i]))  + F2charmM.EvaluatexQ(x[i]/etac, sqrt(q2[i]))  - fthrc*F2charmM0.EvaluatexQ(x[i], sqrt(q2[i]))
	      				       + fthrb*F2bottomZM.EvaluatexQ(x[i], sqrt(q2[i])) + F2bottomM.EvaluatexQ(x[i]/etab, sqrt(q2[i])) - fthrb*F2bottomM0.EvaluatexQ(x[i], sqrt(q2[i]))
	      				       );
	      _flfonll[termID][i] =           (+ FLlight.EvaluatexQ(x[i], sqrt(q2[i]))
	      				       + fthrc*FLcharmZM.EvaluatexQ(x[i], sqrt(q2[i]))  + FLcharmM.EvaluatexQ(x[i]/etac, sqrt(q2[i]))  - fthrc*FLcharmM0.EvaluatexQ(x[i], sqrt(q2[i]))
	      				       + fthrb*FLbottomZM.EvaluatexQ(x[i], sqrt(q2[i])) + FLbottomM.EvaluatexQ(x[i]/etab, sqrt(q2[i])) - fthrb*FLbottomM0.EvaluatexQ(x[i], sqrt(q2[i]))
	      				       );
	      _f3fonll[termID][i] = -charge * (F3light.EvaluatexQ(x[i], sqrt(q2[i])) + F3charmZM.EvaluatexQ(x[i], sqrt(q2[i])) + F3bottomZM.EvaluatexQ(x[i], sqrt(q2[i])));
	      if (!tabulate)
		{
		  _f2fonll[termID][i] =           (+ F2.at(1).Evaluate(x[i], sqrt(q2[i])) + F2.at(2).Evaluate(x[i], sqrt(q2[i])) + F2.at(3).Evaluate(x[i], sqrt(q2[i]))
						   + fthrc*F2.at(4).Evaluate(x[i], sqrt(q2[i]))  + F2M3F.at(4).Evaluate(x[i]/etac, sqrt(q2[i]))  - fthrc*F2M03F.at(4).Evaluate(x[i], sqrt(q2[i]))
						   + fthrb*F2.at(5).Evaluate(x[i], sqrt(q2[i])) + F2M4F.at(5).Evaluate(x[i]/etab, sqrt(q2[i])) - fthrb*F2M04F.at(5).Evaluate(x[i], sqrt(q2[i]))
						   );
		  _flfonll[termID][i] =           (+ FL.at(1).Evaluate(x[i], sqrt(q2[i])) + FL.at(2).Evaluate(x[i], sqrt(q2[i])) + FL.at(3).Evaluate(x[i], sqrt(q2[i]))
						   + fthrc*FL.at(4).Evaluate(x[i], sqrt(q2[i]))  + FLM3F.at(4).Evaluate(x[i]/etac, sqrt(q2[i]))  - fthrc*FLM03F.at(4).Evaluate(x[i], sqrt(q2[i]))
						   + fthrb*FL.at(5).Evaluate(x[i], sqrt(q2[i])) + FLM4F.at(5).Evaluate(x[i]/etab, sqrt(q2[i])) - fthrb*FLM04F.at(5).Evaluate(x[i], sqrt(q2[i]))
						   );
		  _f3fonll[termID][i] = -charge * (F3.at(1).Evaluate(x[i], sqrt(q2[i])) + F3.at(2).Evaluate(x[i], sqrt(q2[i])) + F3.at(3).Evaluate(x[i], sqrt(q2[i])) + F3.at(4).Evaluate(x[i], sqrt(q2[i])) + F3.at(5).Evaluate(x[i], sqrt(q2[i])));
		}

	    }
	}
    }

  bool initcharm = false;
  for (auto termID : _dsIDs)
    if (GetDataFlav(termID) == dataFlav::c)
      initcharm = true;

  if (initcharm)
    {
      const auto F2M  = BuildStructureFunctions(F2ObjM3F,  PDFs3F, PerturbativeOrder, as3F, fBq);
      const auto FLM  = BuildStructureFunctions(FLObjM3F,  PDFs3F, PerturbativeOrder, as3F, fBq);
      const auto F2M0 = BuildStructureFunctions(F2ObjM03F, PDFs3F, PerturbativeOrder, as3F, fBq);
      const auto FLM0 = BuildStructureFunctions(FLObjM03F, PDFs3F, PerturbativeOrder, as3F, fBq);
      const apfel::TabulateObject<apfel::Distribution> F2charmZM {[&] (double const& Q) -> apfel::Distribution{ return F2.at(4).Evaluate(Q);   }, n, qmin, qmax, ord, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> F2charmM  {[&] (double const& Q) -> apfel::Distribution{ return F2M.at(4).Evaluate(Q);  }, n, qmin, qmax, ord, Thresholds3F};
      const apfel::TabulateObject<apfel::Distribution> F2charmM0 {[&] (double const& Q) -> apfel::Distribution{ return F2M0.at(4).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds3F};
      const apfel::TabulateObject<apfel::Distribution> FLcharmZM {[&] (double const& Q) -> apfel::Distribution{ return FL.at(4).Evaluate(Q);   }, n, qmin, qmax, ord, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> FLcharmM  {[&] (double const& Q) -> apfel::Distribution{ return FLM.at(4).Evaluate(Q);  }, n, qmin, qmax, ord, Thresholds3F};
      const apfel::TabulateObject<apfel::Distribution> FLcharmM0 {[&] (double const& Q) -> apfel::Distribution{ return FLM0.at(4).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds3F};
      const apfel::TabulateObject<apfel::Distribution> F3charm {[&] (double const& Q) -> apfel::Distribution{ return F3.at(4).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};

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
	  double mc2 = pow(Thresholds[3],2);

	  for (size_t i = 0; i < Np; i++)
	    {
	      double fthrc = q2[i] > mc2 ? pow(1.-mc2/q2[i],2): 0.;
	      double etac = q2[i] / ( q2[i] + 4 * Masses[3] * Masses[3] );
	      // Skip all points with Q2 below starting scale
	      if (q2[i] < gridLowQLimit*gridLowQLimit)
		continue;
	
	      // Compute structure functions by interpolation in x and Q
	      //	      _f2fonll[termID][i] = F2charm.EvaluatexQ(x[i], sqrt(q2[i]));
	      _f2fonll[termID][i] = fthrc*F2charmZM.EvaluatexQ(x[i], sqrt(q2[i]))  + F2charmM.EvaluatexQ(x[i]/etac, sqrt(q2[i]))  - fthrc*F2charmM0.EvaluatexQ(x[i], sqrt(q2[i]));
	      //              _flfonll[termID][i] = FLcharm.EvaluatexQ(x[i], sqrt(q2[i]));
	      _flfonll[termID][i] = fthrc*FLcharmZM.EvaluatexQ(x[i], sqrt(q2[i]))  + FLcharmM.EvaluatexQ(x[i]/etac, sqrt(q2[i]))  - fthrc*FLcharmM0.EvaluatexQ(x[i], sqrt(q2[i]));
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
      const auto F2M  = BuildStructureFunctions(F2ObjM4F,  PDFs4F, PerturbativeOrder, as4F, fBq);
      const auto FLM  = BuildStructureFunctions(FLObjM4F,  PDFs4F, PerturbativeOrder, as4F, fBq);
      const auto F2M0 = BuildStructureFunctions(F2ObjM04F, PDFs4F, PerturbativeOrder, as4F, fBq);
      const auto FLM0 = BuildStructureFunctions(FLObjM04F, PDFs4F, PerturbativeOrder, as4F, fBq);

      const apfel::TabulateObject<apfel::Distribution> F2bottomZM{[&] (double const& Q) -> apfel::Distribution{ return F2.at(5).Evaluate(Q);   }, n, qmin, qmax, ord, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> F2bottomM {[&] (double const& Q) -> apfel::Distribution{ return F2M.at(5).Evaluate(Q);  }, n, qmin, qmax, ord, Thresholds4F};
      const apfel::TabulateObject<apfel::Distribution> F2bottomM0{[&] (double const& Q) -> apfel::Distribution{ return F2M0.at(5).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds4F};
      const apfel::TabulateObject<apfel::Distribution> FLbottomZM{[&] (double const& Q) -> apfel::Distribution{ return FL.at(5).Evaluate(Q);   }, n, qmin, qmax, ord, Thresholds};
      const apfel::TabulateObject<apfel::Distribution> FLbottomM {[&] (double const& Q) -> apfel::Distribution{ return FLM.at(5).Evaluate(Q);  }, n, qmin, qmax, ord, Thresholds4F};
      const apfel::TabulateObject<apfel::Distribution> FLbottomM0{[&] (double const& Q) -> apfel::Distribution{ return FLM0.at(5).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds4F};
      const apfel::TabulateObject<apfel::Distribution> F3bottom{[&] (double const& Q) -> apfel::Distribution{ return F3.at(5).Evaluate(Q); }, n, qmin, qmax, ord, Thresholds};

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
	  double mb2 = pow(Thresholds[4],2);
	  for (size_t i = 0; i < Np; i++)
	    {
	      // Skip all points with Q2 below starting scale
	      if (q2[i] < gridLowQLimit*gridLowQLimit)
		continue;

	      double fthrb = q2[i] > mb2 ? pow(1.-mb2/q2[i],2): 0.;
	      double etab = q2[i] / ( q2[i] + 4 * Masses[4] * Masses[4] );
	
	      // Compute structure functions by interpolation in x and Q
	      //_f2fonll[termID][i] = F2bottom.EvaluatexQ(x[i], sqrt(q2[i]));
	      //_flfonll[termID][i] = FLbottom.EvaluatexQ(x[i], sqrt(q2[i]));
	      _f2fonll[termID][i] = fthrb*F2bottomZM.EvaluatexQ(x[i], sqrt(q2[i])) + F2bottomM.EvaluatexQ(x[i]/etab, sqrt(q2[i])) - fthrb*F2bottomM0.EvaluatexQ(x[i], sqrt(q2[i]));
	      _flfonll[termID][i] = fthrb*FLbottomZM.EvaluatexQ(x[i], sqrt(q2[i])) + FLbottomM.EvaluatexQ(x[i]/etab, sqrt(q2[i])) - fthrb*FLbottomM0.EvaluatexQ(x[i], sqrt(q2[i]));
	      _f3fonll[termID][i] = -charge * F3bottom.EvaluatexQ(x[i], sqrt(q2[i]));
	    }
	}
    }

  
  
}

// FK mode: rows of F2, FL, xF3 per point, from the evolution operators (FKTables.h). Same
// pieces, Q tabulations and combinations as computeFromPDFs.
void ReactionN3LO_DISNC::buildFK()
{
  using xfitter::fk::OperatorTable;
  using xfitter::fk::TabulateFK;
  const YAML::Node yamlNode = XFITTER_PARS::rootNode["byReaction"]["N3LO_DISNC"];
  const YAML::Node QGrid    = yamlNode["QGrid"];
  bool massive = true;
  if (yamlNode["massive"]) massive = yamlNode["massive"].as<int>();
  int msvordmax = 3;
  if (yamlNode["msvordmax"]) msvordmax = yamlNode["msvordmax"].as<int>();
  const int n = QGrid[0].as<int>();
  const double qmin = QGrid[1].as<double>(), qmax = QGrid[2].as<double>();
  const int ord = QGrid[3].as<int>();
  const int PO = OrderMap(XFITTER_PARS::getParamS("Order")) - 1;
  const int POm = std::min(PO, msvordmax);
  const std::function<std::vector<double>(double const&)> fBq = [] (double const& Q) { return apfel::ElectroWeakCharges(Q, false); };
  const std::function<std::vector<double>(double const&)> fDq = [] (double const& Q) { return apfel::ParityViolatingElectroWeakCharges(Q, false); };
  xfitter::EvolutionAPFELxx* e = _apfelxx;
  const std::function<double(double const&)> as   = [e] (double const& Q) { return e->getAlphaS(Q); };
  const std::function<double(double const&)> as3F = [e] (double const& Q) { return e->getAlphaSNF(3, Q); };
  const std::function<double(double const&)> as4F = [e] (double const& Q) { return e->getAlphaSNF(4, Q); };
  OperatorTable const& E5 = _apfelxx->GetOperatorTable(5);
  OperatorTable const& E3 = _apfelxx->GetOperatorTable(3);
  OperatorTable const& E4 = _apfelxx->GetOperatorTable(4);
  const auto tab5 = [&] (decltype(F2Obj) const& obj, std::vector<int> const& k, decltype(fBq) const& c) {
    return TabulateFK(obj, k, PO, as, c, E5, n, qmin, qmax, ord, Thresholds);
  };
  const auto tab3 = [&] (decltype(F2Obj) const& obj, int order) {
    return TabulateFK(obj, {4}, order, as3F, fBq, E3, n, qmin, qmax, ord, Thresholds3F);
  };
  const auto tab4 = [&] (decltype(F2Obj) const& obj, int order) {
    return TabulateFK(obj, {5}, order, as4F, fBq, E4, n, qmin, qmax, ord, Thresholds4F);
  };
  const int nJoint = static_cast<int>(_apfelxx->InputDistributions().at(0).GetDistributionJointGrid().size());
  _fkRows.clear();

  bool incl = false, charm = false, bottom = false;
  for (auto termID : _dsIDs) {
    incl   = incl   || GetDataFlav(termID) == dataFlav::incl;
    charm  = charm  || GetDataFlav(termID) == dataFlav::c;
    bottom = bottom || GetDataFlav(termID) == dataFlav::b;
  }
  // pieces shared by the dataset types
  std::unique_ptr<OperatorTable> F2light, FLlight, F3light, F2cZM, FLcZM, F3cZM, F2bZM, FLbZM, F3bZM;
  std::unique_ptr<OperatorTable> F2tot, FLtot, F3tot;
  if (incl && !massive) {
    F2tot.reset(new OperatorTable(tab5(F2Obj, {0}, fBq)));
    FLtot.reset(new OperatorTable(tab5(FLObj, {0}, fBq)));
    F3tot.reset(new OperatorTable(tab5(F3Obj, {0}, fDq)));
  }
  if ((incl && massive) || charm || bottom) {
    F2cZM.reset(new OperatorTable(tab5(F2Obj, {4}, fBq)));
    FLcZM.reset(new OperatorTable(tab5(FLObj, {4}, fBq)));
    F3cZM.reset(new OperatorTable(tab5(F3Obj, {4}, fDq)));
    F2bZM.reset(new OperatorTable(tab5(F2Obj, {5}, fBq)));
    FLbZM.reset(new OperatorTable(tab5(FLObj, {5}, fBq)));
    F3bZM.reset(new OperatorTable(tab5(F3Obj, {5}, fDq)));
  }
  if (incl && massive) {
    F2light.reset(new OperatorTable(tab5(F2Obj, {1, 2, 3}, fBq)));
    FLlight.reset(new OperatorTable(tab5(FLObj, {1, 2, 3}, fBq)));
    F3light.reset(new OperatorTable(tab5(F3Obj, {1, 2, 3}, fDq)));
  }

  const double mc2 = pow(Thresholds[3], 2), mb2 = pow(Thresholds[4], 2);
  const auto rowsOf = [&] (unsigned termID) -> FKRows& {
    const size_t Np = GetNpoint(termID);
    FKRows& r = _fkRows[termID];
    r.f2.assign(Np, xfitter::fk::Row(nJoint));
    r.fl.assign(Np, xfitter::fk::Row(nJoint));
    r.f3.assign(Np, xfitter::fk::Row(nJoint));
    r.active.assign(Np, false);
    return r;
  };

  // inclusive datasets
  if (incl) {
    std::unique_ptr<OperatorTable> F2cM, F2cM0, FLcM, FLcM0, F2bM, F2bM0, FLbM, FLbM0;
    if (massive) {
      F2cM.reset(new OperatorTable(tab3(F2ObjM3F, POm)));   F2cM0.reset(new OperatorTable(tab3(F2ObjM03F, POm)));
      FLcM.reset(new OperatorTable(tab3(FLObjM3F, POm)));   FLcM0.reset(new OperatorTable(tab3(FLObjM03F, POm)));
      F2bM.reset(new OperatorTable(tab4(F2ObjM4F, POm)));   F2bM0.reset(new OperatorTable(tab4(F2ObjM04F, POm)));
      FLbM.reset(new OperatorTable(tab4(FLObjM4F, POm)));   FLbM0.reset(new OperatorTable(tab4(FLObjM04F, POm)));
    }
    for (auto termID : _dsIDs) {
      if (GetDataFlav(termID) != dataFlav::incl) continue;
      TermData* td = GetTermData(termID);
      const double gridLowQLimit = td->getPDF()->getQgrid()[0];
      const double charge = GetCharge(termID);
      const auto q2 = *GetBinValues(td, "Q2");
      const auto x = *GetBinValues(td, "x");
      FKRows& r = rowsOf(termID);
      for (size_t i = 0; i < r.active.size(); ++i) {
        if (q2[i] < gridLowQLimit * gridLowQLimit) continue;
        r.active[i] = true;
        const double Q = sqrt(q2[i]);
        if (!massive) {
          r.f2[i].add(*F2tot, x[i], Q, 1);
          r.fl[i].add(*FLtot, x[i], Q, 1);
          r.f3[i].add(*F3tot, x[i], Q, -charge);
          continue;
        }
        const double fthrc = q2[i] > mc2 ? pow(1. - mc2 / q2[i], 2) : 0.;
        const double fthrb = q2[i] > mb2 ? pow(1. - mb2 / q2[i], 2) : 0.;
        const double etac = q2[i] / (q2[i] + 4 * Masses[3] * Masses[3]);
        const double etab = q2[i] / (q2[i] + 4 * Masses[4] * Masses[4]);
        r.f2[i].add(*F2light, x[i], Q, 1);
        r.f2[i].add(*F2cZM, x[i], Q, fthrc);  r.f2[i].add(*F2cM, x[i] / etac, Q, 1);  r.f2[i].add(*F2cM0, x[i], Q, -fthrc);
        r.f2[i].add(*F2bZM, x[i], Q, fthrb);  r.f2[i].add(*F2bM, x[i] / etab, Q, 1);  r.f2[i].add(*F2bM0, x[i], Q, -fthrb);
        r.fl[i].add(*FLlight, x[i], Q, 1);
        r.fl[i].add(*FLcZM, x[i], Q, fthrc);  r.fl[i].add(*FLcM, x[i] / etac, Q, 1);  r.fl[i].add(*FLcM0, x[i], Q, -fthrc);
        r.fl[i].add(*FLbZM, x[i], Q, fthrb);  r.fl[i].add(*FLbM, x[i] / etab, Q, 1);  r.fl[i].add(*FLbM0, x[i], Q, -fthrb);
        r.f3[i].add(*F3light, x[i], Q, -charge);
        r.f3[i].add(*F3cZM, x[i], Q, -charge);
        r.f3[i].add(*F3bZM, x[i], Q, -charge);
      }
    }
  }

  // charm and bottom datasets: massive pieces at the full order, as computeFromPDFs
  for (const int flavour : {4, 5}) {
    const dataFlav type = flavour == 4 ? dataFlav::c : dataFlav::b;
    if (!(flavour == 4 ? charm : bottom)) continue;
    const OperatorTable F2M  = flavour == 4 ? tab3(F2ObjM3F, PO)  : tab4(F2ObjM4F, PO);
    const OperatorTable FLM  = flavour == 4 ? tab3(FLObjM3F, PO)  : tab4(FLObjM4F, PO);
    const OperatorTable F2M0 = flavour == 4 ? tab3(F2ObjM03F, PO) : tab4(F2ObjM04F, PO);
    const OperatorTable FLM0 = flavour == 4 ? tab3(FLObjM03F, PO) : tab4(FLObjM04F, PO);
    OperatorTable const& F2ZM = flavour == 4 ? *F2cZM : *F2bZM;
    OperatorTable const& FLZM = flavour == 4 ? *FLcZM : *FLbZM;
    OperatorTable const& F3ZM = flavour == 4 ? *F3cZM : *F3bZM;
    const double m2 = flavour == 4 ? mc2 : mb2, m = Masses[flavour == 4 ? 3 : 4];
    for (auto termID : _dsIDs) {
      if (GetDataFlav(termID) != type) continue;
      TermData* td = GetTermData(termID);
      const double gridLowQLimit = td->getPDF()->getQgrid()[0];
      const double charge = GetCharge(termID);
      const auto q2 = *GetBinValues(td, "Q2");
      const auto x = *GetBinValues(td, "x");
      FKRows& r = rowsOf(termID);
      for (size_t i = 0; i < r.active.size(); ++i) {
        if (q2[i] < gridLowQLimit * gridLowQLimit) continue;
        r.active[i] = true;
        const double Q = sqrt(q2[i]);
        const double fthr = q2[i] > m2 ? pow(1. - m2 / q2[i], 2) : 0.;
        const double eta = q2[i] / (q2[i] + 4 * m * m);
        r.f2[i].add(F2ZM, x[i], Q, fthr);  r.f2[i].add(F2M, x[i] / eta, Q, 1);  r.f2[i].add(F2M0, x[i], Q, -fthr);
        r.fl[i].add(FLZM, x[i], Q, fthr);  r.fl[i].add(FLM, x[i] / eta, Q, 1);  r.fl[i].add(FLM0, x[i], Q, -fthr);
        r.f3[i].add(F3ZM, x[i], Q, -charge);
      }
    }
  }
}

void ReactionN3LO_DISNC::computeFromFK()
{
  const std::vector<double> f0 = xfitter::fk::Flatten(_apfelxx->InputDistributions());
  for (auto const& t : _fkRows) {
    const FKRows& r = t.second;
    const size_t Np = r.active.size();
    _f2fonll[t.first].resize(Np);
    _flfonll[t.first].resize(Np);
    _f3fonll[t.first].resize(Np);
    for (size_t i = 0; i < Np; ++i) {
      if (!r.active[i]) {
        _f2fonll[t.first][i] = _flfonll[t.first][i] = _f3fonll[t.first][i] = 0;
        continue;
      }
      _f2fonll[t.first][i] = r.f2[i].contract(f0);
      _flfonll[t.first][i] = r.fl[i].contract(f0);
      _f3fonll[t.first][i] = r.f3[i].contract(f0);
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
