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

namespace xfitter
{

// the class factories
extern "C" EvolutionLHAPDF*create(const char*name){
  return new EvolutionLHAPDF(name);
}
EvolutionLHAPDF::EvolutionLHAPDF(const char*name):BaseEvolution(name){
  _pdf=nullptr;
}
const char*EvolutionLHAPDF::getClassName()const{return "LHAPDF";}


/// Global initialization
  void EvolutionLHAPDF::atStart() {
    atConfigurationChange();
  };
  void EvolutionLHAPDF::atConfigurationChange() {
    using namespace std;
    YAML::Node pars=XFITTER_PARS::getEvolutionNode(_name);
    try{
      _set_name =pars["set"].as<std::string>();
    }catch(YAML::TypedBadConversion<std::string>&ex){
      if(!pars["set"]){
        cerr<<"[ERROR] No set name given for LHAPDF Evolution \""<<_name<<"\""<<endl;
        hf_errlog(2018101230,"F: No set name given for LHAPDF Evolution, see stderr");
      }
      cerr<<"[ERROR] Failed to parse set name for LHAPDF Evolution \""<<_name<<"\""<<endl;
      hf_errlog(2018101231,"F: Failed to parse set name for LHAPDF Evolution, see stderr");
    }
    try{
      _member   =pars["member"].as<int>();
    }catch(YAML::TypedBadConversion<int>&ex){
      if(!pars["member"]){
        cerr<<"[ERROR] No member id given for LHAPDF Evolution \""<<_name<<"\""<<endl;
        hf_errlog(2018101232,"F: No member id given for LHAPDF Evolution, see stderr");
      }
      cerr<<"[ERROR] Failed to parse member id for LHAPDF Evolution \""<<_name<<"\""<<endl;
      hf_errlog(2018101233,"F: Failed to parse member id for LHAPDF Evolution, see stderr");
    }
    CheckForPDF(_set_name.c_str());
    if(_pdf)delete _pdf;
    _pdf=LHAPDF::mkPDF(_set_name,_member);
  };

  /// Return PDFs as a map <int,double> where int is PDF ID (-6, ... 6, 21)
  std::function<std::map<int,double>(double const& x, double const& Q)> EvolutionLHAPDF::xfxQMap() {
    auto f0 = [=] (double const& x, double const& Q)->std::map<int, double>{ return _pdf->xfxQ(x, Q); };
    return  f0;
  };

  /// Returns PDFs as a function of i, x, Q
  std::function<double(int const& i, double const& x, double const& Q)> EvolutionLHAPDF::xfxQDouble() {
    auto f0 = [=] (int const& i, double const& x, double const& Q)->double { return _pdf->xfxQ(i,x, Q); };
    return  f0;
  };

  /// Returns PDFs as double pdfs* --> double[13] from -6 to 6.
  std::function<void(double const& x, double const& Q, double* pdfs)> EvolutionLHAPDF::xfxQArray() {
    auto f0 = [=] (double const& x, double const& Q, double* pdfs)->void
      {
        std::vector<double> vpdfs{13};
        _pdf->xfxQ(x, Q, vpdfs);
        for (int i=0; i<13; i++)
          pdfs[i] = vpdfs[i];
        return ;
      };
    return  f0;
  };

  /// Returns alphaS
  std::function<double(double const& Q)> EvolutionLHAPDF::AlphaQCD() {
    auto f0 = [=] ( double const& Q) -> double {
      return _pdf->alphasQ(Q);
    };
    return f0;
  };

  /// Get property
  std::string EvolutionLHAPDF::getPropertyS(std::string const& propertyName) const {
    if ( _pdf->info().has_key(propertyName) ) {
      return _pdf->info().get_entry(propertyName);
    }
    else {
      hf_errlog(2018082421,"S: Missing PDF information for property "+propertyName);
      return "";
    };
  }

  /// Get property
  int EvolutionLHAPDF::getPropertyI(std::string const& propertyName) const {
    std::string sVal = getPropertyS(propertyName);
    return std::stoi(sVal);
  }

}
