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
    _mch = XFITTER_PARS::getParamD("mch");
    _mbt = XFITTER_PARS::getParamD("mbt");
    _mtp = XFITTER_PARS::getParamD("mtp");
    _mch_last = *_mch;
    _mbt_last = *_mbt;
    _mtp_last = *_mtp;
    const YAML::Node xGrid = yamlNode["xGrid"];

    vector<apfel::SubGrid> sgv;
    for(auto const& sg : xGrid)
      sgv.push_back(apfel::SubGrid{sg[0].as<int>(), sg[1].as<double>(), sg[2].as<int>()});

    // x-space grid (the grid parameters should be in parameters.yaml
    // in APFELxx/yaml. I need to find a clever way to retrieve them)
    _Grid = std::unique_ptr<const apfel::Grid>(new apfel::Grid(sgv));

    // Vectors of masses and thresholds
    _Masses = {0, 0, 0, *_mch, *_mbt, *_mtp};
    _Thresholds = _Masses;
    _isFFNS = 0; // VFNS by default
    if(XFITTER_PARS::gParametersI.find("isFFNS") != XFITTER_PARS::gParametersI.end()) {
      _isFFNS = XFITTER_PARS::gParametersI.at("isFFNS");
    }
    if (yamlNode["isFFNS"]) {
      _isFFNS = yamlNode["isFFNS"].as<int>();
    }
    _NFlavour = -1;
    _NFlavour = XFITTER_PARS::gParametersI.at("NFlavour");
    if (yamlNode["NFlavour"]) {
      _NFlavour = yamlNode["NFlavour"].as<int>();
    }
    if(_isFFNS == 1) {
      for (int i = 5; i >= _NFlavour; i--) {
        _Thresholds[i] = 1.0e10;
      }
    }
    else if(_isFFNS == 0) {;}
    else if(_isFFNS != 0) {
      hf_errlog(2025020101, "F: Unsupported _isFFNS = " + std::to_string(_isFFNS));
    }
    if(XFITTER_PARS::gParametersS.find("heavyQuarkMassScheme") != XFITTER_PARS::gParametersS.end()) {
      _heavyQuarkMassScheme = XFITTER_PARS::gParametersS.at("heavyQuarkMassScheme");
    }
    if (yamlNode["heavyQuarkMassScheme"]) {
      _heavyQuarkMassScheme = yamlNode["heavyQuarkMassScheme"].as<string>();
    }
    // Initialize QCD evolution objects
    if (_heavyQuarkMassScheme == "Pole") {
      _DglapObj = apfel::InitializeDglapObjectsQCD(*_Grid, _Masses, _Thresholds);
    }
    else if (_heavyQuarkMassScheme == "MSBar") {
      _DglapObj = apfel::InitializeDglapObjectsQCDMSbarMass(*_Grid, _Masses, _Thresholds);
    }
    else {
      hf_errlog(2025020501, "F: Unsupported _heavyQuarkMassScheme = " + _heavyQuarkMassScheme);
    }
    //_DglapObj.FlavourScheme = FFNS;
    atConfigurationChange();
  }

  //_________________________________________________________________________________
  void EvolutionAPFELxx::atIteration()
  {
    const YAML::Node yamlNode=XFITTER_PARS::getEvolutionNode(_name);
    // Restart from scratch at each iteration, e.g. when fitting heavy-quark masses (fast enough)
    if (yamlNode["restart_at_each_iteration"] && yamlNode["restart_at_each_iteration"].as<int>() == 1) {
      // Restart only if any of relevant heavy quark masses changed
      if (_isFFNS == 0) {
        if (_NFlavour >= 4 && _mch_last != *_mch || _NFlavour >= 5 && _mbt_last != *_mbt || _NFlavour >= 6 &&  _mtp_last != *_mtp) {
          this->atStart();
        }
      }
    }
    // Retrieve the relevant parameters needed to compute the evolutions
    const int     PtOrder    = OrderMap(XFITTER_PARS::getParamS("Order")) - 1;
    _Q0 = *XFITTER_PARS::getParamD((yamlNode["Q0"]) ? yamlNode["Q0"].as<string>() : "Q0");
    const double *Mz = XFITTER_PARS::getParamD("Mz");
    try {
      _alphas_q0 = XFITTER_PARS::getParamD((yamlNode["alphas_Q0"]) ? yamlNode["alphas_Q0"].as<string>() : "alphas_Q0");
    }
    catch(std::out_of_range&ex) {
      _alphas_q0 = nullptr;
    }
    if (!_alphas_q0) {
      _alphas_q0 = Mz;
    }
    _alphas = XFITTER_PARS::getParamD((yamlNode["alphas"]) ? yamlNode["alphas"].as<string>() : "alphas");
    const YAML::Node QGrid   = yamlNode["QGrid"];

    // Reinitialise and tabulate the running coupling at every
    // iteration. This is fast enough and allows for the reference
    // value to be fitted.
    if (_heavyQuarkMassScheme == "Pole") {
      apfel::AlphaQCD a{*_alphas, *_alphas_q0, _Masses, _Thresholds, PtOrder};
      const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
      _AlphaQCD = [=] (double const& mu) -> double{ return Alphas.Evaluate(mu); };
    }
    else if (_heavyQuarkMassScheme == "MSBar") {
      apfel::AlphaQCDMSbarMass a{*_alphas, *_alphas_q0, _Masses, _Thresholds, PtOrder};
      const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
      _AlphaQCD = [=] (double const& mu) -> double{ return Alphas.Evaluate(mu); };
    }

    // Construct the DGLAP objects
    const auto Dglap = BuildDglap(_DglapObj,
      [=] (double const& x, double const&)->std::map<int,double>{
        return apfel::PhysToQCDEv(_inPDFs->xfxMap(x));
      },
    _Q0, PtOrder, _AlphaQCD);

    // Tabulate PDFs (ideally the parameters of the tabulation should
    // be read from parameters.yaml).
    _TabulatedPDFs = std::unique_ptr<apfel::TabulateObject<apfel::Set<apfel::Distribution>>>
      (new apfel::TabulateObject<apfel::Set<apfel::Distribution>>{*Dglap,
          QGrid[0].as<int>(),
          QGrid[1].as<double>(),
          QGrid[2].as<double>(),
          QGrid[3].as<int>()});
  }
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
