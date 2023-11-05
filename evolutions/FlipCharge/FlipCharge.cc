#include"FlipCharge.h"
#include"xfitter_cpp_base.h"
#include"xfitter_pars.h"
#include"xfitter_steer.h"
#include"iostream"
#include<utility> //for swap
#include <map>

namespace xfitter {
  extern"C" FlipChargeEvol* create(const char* s) {
    return new FlipChargeEvol(s);
  } //for dynamic loading

  void FlipChargeEvol::atStart() {
    atConfigurationChange();
  }

  void FlipChargeEvol::atConfigurationChange() {
    const YAML::Node node=XFITTER_PARS::getEvolutionNode(_name)["input"];
    try {
      input = get_evolution(node.as<std::string>());
    } catch (const YAML::BadConversion& ex) {
      if(!node) {
        //If no input is provided, use default evolution
        input = defaultEvolutionInstance();
        std::cerr << "[INFO] Evolution \"" << _name << "\" of class " << getClassName()
                  << ": \"input\" not provided; using default evolution \"" << input->_name << "\"" << std::endl;
      } else {
        std::cerr << "[ERROR] In evolution " << _name << " of class "<< getClassName()
                  << ": failed to convert input to string, YAML node:" << node << std::endl;
        hf_errlog(18122501, "F: Bad input to evolution FlipCharge, see stderr");
      }
    }
  }

  std::map<int, double>FlipChargeEvol::xfxQmap(double x,double Q) {
    double p[13];
    input->xfxQarray(x,Q,p);
    return {
      {-6, p[12]},
      {-5, p[11]},
      {-4, p[10]},
      {-3, p[ 9]},
      {-2, p[ 8]},
      {-1, p[ 7]},
      {21, p[ 6]},
      { 1, p[ 5]},
      { 2, p[ 4]},
      { 3, p[ 3]},
      { 4, p[ 2]},
      { 5, p[ 1]},
      { 6, p[ 0]}
    };
  }

  void FlipChargeEvol::xfxQarray(double x, double Q, double* p) {
    input->xfxQarray(x,Q,p);
    std::swap(p[0],p[12]);
    std::swap(p[1],p[11]);
    std::swap(p[2],p[10]);
    std::swap(p[3],p[ 9]);
    std::swap(p[4],p[ 8]);
    std::swap(p[5],p[ 7]);
  }

  double FlipChargeEvol::xfxQ(int i,double x,double Q) {
    return input->xfxQ(-i, x, Q);
  }

  double FlipChargeEvol::getAlphaS(double Q) {
    return input->getAlphaS(Q);
  }

  std::vector<double> FlipChargeEvol::getXgrid() {
    return input->getXgrid();
  }

  std::vector<double> FlipChargeEvol::getQgrid() {
    return input->getQgrid();
  }
}
