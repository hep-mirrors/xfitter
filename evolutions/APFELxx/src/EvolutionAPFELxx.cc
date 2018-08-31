#include "EvolutionAPFELxx.h"
#include "xfitter_pars.h"
#include "xfitter_cpp_base.h"

#include <apfel/alphaqcd.h>
#include <apfel/messages.h>
#include <apfel/rotations.h>


namespace xfitter
{
  // the class factories
  extern "C" EvolutionAPFELxx*create(const char*name){
    return new EvolutionAPFELxx(name);
  }


  //_________________________________________________________________________________
  void EvolutionAPFELxx::initFromYaml(const YAML::Node yamlNode)
  {
    // APFEL++ banner
    apfel::Banner();

		_inPDFs=XFITTER_PARS::getInputFunctionFromYaml(yamlNode);
    // Retrieve parameters needed to initialize APFEL++.
    const double* MCharm   = XFITTER_PARS::gParameters.at("mch");
    const double* MBottom  = XFITTER_PARS::gParameters.at("mbt");
    const double* MTop     = XFITTER_PARS::gParameters.at("mtp");
    const YAML::Node xGrid = yamlNode["xGrid"];

    vector<apfel::SubGrid> sgv;
    for(auto const& sg : xGrid)
      sgv.push_back(apfel::SubGrid{sg[0].as<int>(), sg[1].as<double>(), sg[2].as<int>()});

    // x-space grid (the grid parameters should be in parameters.yaml
    // in APFELxx/yaml. I need to find a clever way to retrieve them)
    _Grid = std::unique_ptr<const apfel::Grid>(new apfel::Grid(sgv));

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
    const int     PtOrder    = OrderMap(XFITTER_PARS::gParametersS.at("Order")) - 1;
    const double* Q0         = XFITTER_PARS::gParameters.at("Q0");
    const double* Q_ref      = XFITTER_PARS::gParameters.at("Mz");
    const double* Alphas_ref = XFITTER_PARS::gParameters.at("alphas");
		//XXX HACKS XXX
		//I do not understand how QGrid is supposed to change between iterations --Ivan
		//This will not work with the new parameters.yaml syntax
		//TODO
    const YAML::Node QGrid   = XFITTER_PARS::gParametersY.at("APFELxx")["QGrid"];

    // Reinitialise and tabulate the running coupling at every
    // iteration. This is fast enough and allows for the reference
    // value to be fitted.
    apfel::AlphaQCD a{*Alphas_ref, *Q_ref, _Masses, _Thresholds, PtOrder};
    const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
    _AlphaQCD = [=] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

    // Construct the DGLAP objects
    const auto Dglap = BuildDglap(_DglapObj,
				  [=] (double const& x, double const&)->std::map<int,double>{ return apfel::PhysToQCDEv(this->_inPDFs(x)); },
				  *Q0, PtOrder, _AlphaQCD);

    // Tabulate PDFs (ideally the parameters of the tabulation should
    // be read from parameters.yaml).
    _TabulatedPDFs = std::unique_ptr<apfel::TabulateObject<apfel::Set<apfel::Distribution>>>
      (new apfel::TabulateObject<apfel::Set<apfel::Distribution>>{*Dglap,
	  QGrid[0].as<int>(),
	  QGrid[1].as<double>(),
	  QGrid[2].as<double>(),
	  QGrid[3].as<int>()});
  }

  //_________________________________________________________________________________
  std::function<std::map<int,double>(double const& x, double const& Q)> EvolutionAPFELxx::xfxQMap()
  {
    // return lambda function straight away.
    return [=] (double const& x, double const& Q) -> std::map<int,double>{ return _TabulatedPDFs->EvaluateMapxQ(x, Q); };
  }

  //_________________________________________________________________________________
  std::function<double(int const& i, double const& x, double const& Q)> EvolutionAPFELxx::xfxQDouble()
  {
    // return lambda function straight away.
    return [=] (int const& i, double const& x, double const& Q) -> double{ return _TabulatedPDFs->EvaluatexQ(i, x, Q); };
  }

  //_________________________________________________________________________________
  std::function<void(double const& x, double const& Q, double* pdfs)> EvolutionAPFELxx::xfxQArray()
  {
    // return lambda function straight away.
    return [=] (double const& x, double const& Q, double* pdfs) -> void
      {
	// Get map of PDFs
	const std::map<int,double> fset = apfel::QCDEvToPhys(_TabulatedPDFs->EvaluateMapxQ(x, Q));

	// Fill in array of PDFs to be returned
	//       	int counter = 0;
	//for(auto const& f : fset) 
	//  pdfs[counter++] = f.second;

	pdfs[0] = fset.at(-6);
	pdfs[1] = fset.at(-5);
	pdfs[2] = fset.at(-4);
	pdfs[3] = fset.at(-3);
	pdfs[4] = fset.at(-2);
	pdfs[5] = fset.at(-1);
	pdfs[6] = fset.at(0);
	pdfs[7] = fset.at(1);
	pdfs[8] = fset.at(2);
	pdfs[9] = fset.at(3);
	pdfs[10] = fset.at(4);
	pdfs[11] = fset.at(5);
	pdfs[12] = fset.at(6);
      };
  }
}
