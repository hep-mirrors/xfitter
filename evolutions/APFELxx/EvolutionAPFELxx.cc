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
  const char*EvolutionAPFELxx::getClassName()const{return "APFELxx";}


  //_________________________________________________________________________________
  void EvolutionAPFELxx::atStart()
  {
    // APFEL++ banner
    apfel::Banner();

    const YAML::Node yamlNode=XFITTER_PARS::getEvolutionNode(_name);
    _inPDFs=XFITTER_PARS::getInputDecomposition(yamlNode);
    // Retrieve parameters needed to initialize APFEL++.
    const double* MCharm   = XFITTER_PARS::getParamD("mch");
    const double* MBottom  = XFITTER_PARS::getParamD("mbt");
    const double* MTop     = XFITTER_PARS::getParamD("mtp");
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
    //_DglapObj = apfel::InitializeDglapObjectsQCD(*_Grid, _Masses, _Thresholds);
    std::vector<int>  IMod = {0, 0, 0, 0, 0, 0, 0};
    //imod 1 : A
    //imod 2 : B
    //imod 0 : (A+B)/2
    IMod[0] = yamlNode["P3NSp"].as<int>();
    IMod[1] = yamlNode["P3NSm"].as<int>();
    IMod[2] = yamlNode["P3NSs"].as<int>();
    IMod[3] = yamlNode["P3SGps"].as<int>();
    IMod[4] = yamlNode["P3SGqg"].as<int>();
    IMod[5] = yamlNode["P3SGgq"].as<int>();
    IMod[6] = yamlNode["P3SGgg"].as<int>();    
    _DglapObj = apfel::InitializeDglapObjectsQCD(*_Grid, _Masses, _Thresholds, false, 1e-5, IMod);
    atConfigurationChange();
  }

  //_________________________________________________________________________________
  void EvolutionAPFELxx::atIteration()
  {
    const YAML::Node yamlNode=XFITTER_PARS::getEvolutionNode(_name);
    // Retrieve the relevant parameters needed to compute the evolutions
    const int     PtOrder    = OrderMap(XFITTER_PARS::getParamS("Order")) - 1;
    const double* Q0         = XFITTER_PARS::getParamD("Q0");
    const double* Q_ref      = XFITTER_PARS::getParamD("Mz");
    const double* Alphas_ref = XFITTER_PARS::getParamD("alphas");
    const YAML::Node QGrid   = yamlNode["QGrid"];

    // Reinitialise and tabulate the running coupling at every
    // iteration. This is fast enough and allows for the reference
    // value to be fitted.
    apfel::AlphaQCD a{*Alphas_ref, *Q_ref, _Masses, _Thresholds, PtOrder};
    const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
    _AlphaQCD = [=] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

    // Construct the DGLAP objects
    const auto Dglap = BuildDglap(_DglapObj,
      [=] (double const& x, double const&)->std::map<int,double>{
        return apfel::PhysToQCDEv(_inPDFs->xfxMap(x));
      },
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
  apfel::TabulateObject<apfel::Set<apfel::Distribution>> EvolutionAPFELxx::GetTabulatedPDFs() {return *_TabulatedPDFs;};
  std::function<double(double const& Q)>                 EvolutionAPFELxx::GetAlphaQCD() {return _AlphaQCD;};
  
  std::map<int,double>EvolutionAPFELxx::xfxQmap(double x,double Q){
    return apfel::QCDEvToPhys(_TabulatedPDFs->EvaluateMapxQ(x,Q));
  }
  double EvolutionAPFELxx::xfxQ(int i,double x,double Q){
    return _TabulatedPDFs->EvaluatexQ(i,x,Q);
  }
  void EvolutionAPFELxx::xfxQarray(double x,double Q,double*pdfs){
    // Get map of PDFs
    const std::map<int,double> fset = apfel::QCDEvToPhys(_TabulatedPDFs->EvaluateMapxQ(x, Q));
    pdfs[0] =fset.at(-6);
    pdfs[1] =fset.at(-5);
    pdfs[2] =fset.at(-4);
    pdfs[3] =fset.at(-3);
    pdfs[4] =fset.at(-2);
    pdfs[5] =fset.at(-1);
    pdfs[6] =fset.at(0);
    pdfs[7] =fset.at(1);
    pdfs[8] =fset.at(2);
    pdfs[9] =fset.at(3);
    pdfs[10]=fset.at(4);
    pdfs[11]=fset.at(5);
    pdfs[12]=fset.at(6);
  }
  double EvolutionAPFELxx::getAlphaS(double Q){
    return _AlphaQCD(Q);
  }

  vector<double> EvolutionAPFELxx::getXgrid() {
    return _Grid->GetJointGrid().GetGrid();
  }

  vector<double> EvolutionAPFELxx::getQgrid() {
    return _TabulatedPDFs->GetQGrid();
  }
}
