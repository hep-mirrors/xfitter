#include "EvolutionAPFELxx.h"
#include "xfitter_pars.h"

#include <apfel/alphaqcd.h>
#include <apfel/messages.h>
#include <apfel/lhtoypdfs.h>

namespace xfitter
{
  //_________________________________________________________________________________
  void EvolutionAPFELxx::initAtStart()
  {
    // APFEL++ banner
    apfel::Banner();

    // Retrieve parameters needed to initialize APFEL++.
    const double* MCharm     = XFITTER_PARS::gParameters.at("mch");
    const double* MBottom    = XFITTER_PARS::gParameters.at("mbt");
    const double* MTop       = XFITTER_PARS::gParameters.at("mtp");

    // x-space grid (the grid parameters should be in parameters.yaml
    // in APFELxx/yaml. I need to find a clever way to retrieve them)
    _Grid = std::unique_ptr<const apfel::Grid>(new apfel::Grid({
	  apfel::SubGrid{100, 1e-5, 3},
	    apfel::SubGrid{60, 1e-1, 3},
	      apfel::SubGrid{50, 6e-1, 3},
		apfel::SubGrid{50, 8e-1, 3}}));

    // Vectors of masses and thresholds
    _Masses = {0, 0, 0, *MCharm, *MBottom, *MTop};
    _Thresholds = _Masses;

    // Initialize QCD evolution objects
    _DglapObj = apfel::InitializeDglapObjectsQCD(*_Grid, _Masses, _Thresholds);
  }

  //_________________________________________________________________________________
  void EvolutionAPFELxx::initAtIteration()
  {
    // Retrieve the relevant parameters needed to compute the evolutions
    const int     PtOrder    = XFITTER_PARS::gParametersI.at("Order") - 1;
    const double* Q_ref      = XFITTER_PARS::gParameters.at("Mz");
    const double* Alphas_ref = XFITTER_PARS::gParameters.at("alphas");

    // Reinitialise and tabulate the running coupling at every
    // iteration. This is fast enough and allows for the reference
    // value to be fitted.
    apfel::AlphaQCD a{*Alphas_ref, *Q_ref, _Masses, _Thresholds, PtOrder};
    const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
    const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

    // Construct the DGLAP objects
    _Dglap = BuildDglap(_DglapObj, apfel::LHToyPDFs, 1, PtOrder, as);

    // Tabulate PDFs (ideally the parameters of the tabulation should
    // be read from parameters.yaml).
    _TabulatedPDFs = std::unique_ptr<apfel::TabulateObject<apfel::Set<apfel::Distribution>>>
      (new apfel::TabulateObject<apfel::Set<apfel::Distribution>>{*_Dglap, 50, 1, 1000, 3});
  }
}
