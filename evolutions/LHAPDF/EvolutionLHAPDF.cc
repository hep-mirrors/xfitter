/*
   @file EvolutionLHAPDF.cc
   @date  2018-08-20
   @author  AddEvolution.py
   Created by  AddEvolution.py on 2018-08-20
*/

#include "EvolutionLHAPDF.h"
#include "xfitter_pars.h"
#include "CheckForPDF.h"
#include "xfitter_cpp_base.h"

namespace xfitter {

// the class factories
  extern "C" EvolutionLHAPDF* create(const char* name)
  {
    return new EvolutionLHAPDF(name);
  }
  EvolutionLHAPDF::EvolutionLHAPDF(const char* name): BaseEvolution(name)
  {
    _pdf = nullptr;
  }
  const char* EvolutionLHAPDF::getClassName()const {return "LHAPDF";}


/// Global initialization
  void EvolutionLHAPDF::atStart()
  {
    LHAPDF::Info& cfg = LHAPDF::getConfig();
    cfg.set_entry("Verbosity", 0);
    atConfigurationChange();
  };

  void EvolutionLHAPDF::atConfigurationChange()
  {
    using namespace std;
    YAML::Node pars = XFITTER_PARS::getEvolutionNode(_name);
    try {
      _set_name = pars["set"].as<std::string>();
    } catch (YAML::TypedBadConversion<std::string>& ex) {
      if (!pars["set"]) {
        cerr << "[ERROR] No set name given for LHAPDF Evolution \"" << _name << "\"" << endl;
        hf_errlog(2018101230, "F: No set name given for LHAPDF Evolution, see stderr");
      }
      cerr << "[ERROR] Failed to parse set name for LHAPDF Evolution \"" << _name << "\"" << endl;
      hf_errlog(2018101231, "F: Failed to parse set name for LHAPDF Evolution, see stderr");
    }
    try {
      _member   = pars["member"].as<int>();
    } catch (YAML::TypedBadConversion<int>& ex) {
      if (!pars["member"]) {
        cerr << "[ERROR] No member id given for LHAPDF Evolution \"" << _name << "\"" << endl;
        hf_errlog(2018101232, "F: No member id given for LHAPDF Evolution, see stderr");
      }
      cerr << "[ERROR] Failed to parse member id for LHAPDF Evolution \"" << _name << "\"" << endl;
      hf_errlog(2018101233, "F: Failed to parse member id for LHAPDF Evolution, see stderr");
    }
    CheckForPDF(_set_name.c_str());
    if (_pdf)delete _pdf;
    _pdf = LHAPDF::mkPDF(_set_name, _member);

    // Also get some extra info:
    auto const& pdfs = _pdf->xfxQ(0.1, 10.0);
    _has_photon = pdfs.count(22) > 0;
  };

  /// Return PDFs as a map <int,double> where int is PDF ID (-6, ... 6, 21)
  std::map<int, double>EvolutionLHAPDF::xfxQmap(double x, double Q)
  {
    return _pdf->xfxQ(x, Q);
  };

  /// Returns PDFs as a function of i, x, Q
  double EvolutionLHAPDF::xfxQ(int i, double x, double Q)
  {
    return _pdf->xfxQ(i, x, Q);
  };

  /// Returns PDFs as double pdfs* --> double[13] from -6 to 6.
  void EvolutionLHAPDF::xfxQarray(double x, double Q, double* pdfs)
  {
    // some APPLgrids have values x = 0 or x = 1, LHAPDF throws LHAPDF::RangeError
    if (x <= 0.0 || x >= 1.0) {
      for (int i = 0; i < 13; i++)
        pdfs[i] = 0.0;
      return;
    }
    std::vector<double> vpdfs{13};
    _pdf->xfxQ(x, Q, vpdfs);
    for (int i = 0; i < 13; i++)
      pdfs[i] = vpdfs[i];

    if (_has_photon) {
      pdfs[13] = _pdf->xfxQ(22, x, Q);
    }

    return ;
  };

  /// Returns alphaS
  double EvolutionLHAPDF::getAlphaS(double Q)
  {
    return _pdf->alphasQ(Q);
  };

  /// Get property
  std::string EvolutionLHAPDF::getPropertyS(std::string const& propertyName) const
  {
    if (_pdf->info().has_key(propertyName)) {
      return _pdf->info().get_entry(propertyName);
    } else {
      hf_errlog(2018082421, "S: Missing PDF information for property " + propertyName);
      return "";
    };
  }

  /// Get property
  int EvolutionLHAPDF::getPropertyI(std::string const& propertyName) const
  {
    std::string sVal = getPropertyS(propertyName);
    return std::stoi(sVal);
  }

  /// Get property
  double EvolutionLHAPDF::getPropertyD(std::string const& propertyName) const {
    std::string sVal = getPropertyS(propertyName);
    return std::stod(sVal);
  }
  
  /// Get property with default value
  double EvolutionLHAPDF::getPropertyD(std::string const& propertyName, double defval) const {
      double val = _pdf->info().get_entry_as<double>(propertyName, -1);
      return val;
  }
}
