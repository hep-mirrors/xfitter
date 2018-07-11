#include "EvolutionAPFELxx.h"
#include "xfitter_pars.h"

#include <apfel/grid.h>
#include <apfel/alphaqcd.h>
#include <apfel/dglap.h>
#include <apfel/dglapbuilder.h>
#include <apfel/tabulateobject.h>

namespace xfitter
{
  void EvolutionAPFELxx::initAtStart() const
  {
    // Retrieve parameters needed to initialize APFEL++.
    const int     PtOrder    = XFITTER_PARS::gParametersI.at("Order") - 1; //OrderMap(GetParamS("Order")) - 1;
    const double* MCharm     = XFITTER_PARS::gParameters.at("mch"); //GetParam("mch");
    const double* MBottom    = XFITTER_PARS::gParameters.at("mbt"); //GetParam("mbt");
    const double* MTop       = XFITTER_PARS::gParameters.at("mtp"); //GetParam("mtp");
    const double* Q_ref      = XFITTER_PARS::gParameters.at("Mz"); //GetParam("Mz");;
    const double* Alphas_ref = XFITTER_PARS::gParameters.at("alphas"); //GetParam("alphas");

    // x-space grid
    const apfel::Grid g{{apfel::SubGrid{100,1e-5,3}, apfel::SubGrid{60,1e-1,3}, apfel::SubGrid{50,6e-1,3}, apfel::SubGrid{50,8e-1,3}}};

    // Vectors of masses and thresholds
    const vector<double> Masses = {0, 0, 0, *MCharm, *MBottom, *MTop};
    const vector<double> Thresholds = Masses;

    // Initialize QCD evolution objects
    const auto DglapObj = apfel::InitializeDglapObjectsQCD(g, Masses, Thresholds);

    // Running coupling
    apfel::AlphaQCD a{*Alphas_ref, *Q_ref, Masses, PtOrder};
    const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
    const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

    // Construct the DGLAP objects
    //auto EvolvedPDFs = BuildDglap(DglapObj, LHToyPDFs, mu0, PerturbativeOrder, as);
  }
}
