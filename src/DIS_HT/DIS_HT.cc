#include "DIS_HT.h"
#include "spline.h"
#include "TermData.h"
#include <cmath>
#include <sstream>
#include <iostream>

DIS_HT::DIS_HT(TermData* td) {
  _spline_f2 = new tk::spline();
  _spline_ft = new tk::spline();
  // for arrays, expect comma-separated strings
  // each item should be either double or parameter name
  auto read_array = [td](const std::string& parname, std::vector<bool>& isnew) {
    isnew.clear();
    std::istringstream ss(td->getParamS(parname));
    std::string token;
    int counter = 0;
    std::vector<const double* > result;
    while(std::getline(ss, token, ','))
    {
      if (td->hasParam(token)) {
        std::cout << "DIS_HT: using higher twist parameter " + parname + "[" + std::to_string(counter) + "] = \"" + token + "\" as xfitter parameter" << std::endl;
        result.push_back(td->getParamD(token));
        isnew.push_back(false);
      }
      else {
        try {
          std::cout << "DIS_HT: using higher twist parameter " + parname + "[" + std::to_string(counter) + "] = \"" + token + "\" as constant value" << std::endl;
          result.push_back(new double(stod(token)));
          isnew.push_back(true);
        }
        catch (const std::invalid_argument&) {
          string msg = "F: DIS_HT requested higher twist parameter " + parname + "[" + std::to_string(counter) + "] = \"" + token + "\" is not a parameter name or double";
          hf_errlog(24051700, msg);
        }
      }
      counter++;
    }
    //printf("read_array result ="); for(auto& r : result) printf(" %.3f", *r); printf("\n");
    return result;
  };
  auto read_double = [td](const std::string& parname, bool& isnew) {
    std::istringstream ss(td->getParamS(parname));
    std::string token;
    const double* result;
    while(std::getline(ss, token, ','))
    {
      if (td->hasParam(token)) {
        std::cout << "DIS_HT: using higher twist parameter " + parname + " = \"" + token + "\" as xfitter parameter" << std::endl;
        result = td->getParamD(token);
        isnew = false;
      }
      else {
        try {
          std::cout << "DIS_HT: using higher twist parameter " + parname + " = \"" + token + "\" as constant value" << std::endl;
          result = new double(stod(token));
          isnew = true;
        }
        catch (const std::invalid_argument&) {
          string msg = "F: DIS_HT requested higher twist parameter " + parname + " = \"" + token + "\" is not a parameter name or double";
          hf_errlog(24051700, msg);
        }
      }
    }
    //printf("read_double result = %.3f\n", *result);
    return result;
  };
  _ht_x = read_array("ht_x", _isnew_ht_x);
  _ht_2 = read_array("ht_vals_f2", _isnew_ht_2);
  _ht_t = read_array("ht_vals_ft", _isnew_ht_t);
  _ht_alpha_2 = read_double("ht_val_alpha_f2", _isnew_ht_alpha_2);
  _ht_alpha_t = read_double("ht_val_alpha_ft", _isnew_ht_alpha_t);
  if (td->hasParam("ht_val_scale")) {
    _ht_scale = read_double("ht_val_scale", _isnew_ht_scale);
  }
  else {
    _ht_scale = new double(1);
    _isnew_ht_scale = true;
  }
}

DIS_HT::~DIS_HT() {
  for(size_t i = 0; i < _isnew_ht_x.size(); i++) {
    if(_isnew_ht_x[i]) {
      delete _ht_x[i];
    }
  }
  for(size_t i = 0; i < _isnew_ht_2.size(); i++) {
    if(_isnew_ht_2[i]) {
      delete _ht_2[i];
    }
  }
  for(size_t i = 0; i < _isnew_ht_t.size(); i++) {
    if(_isnew_ht_t[i]) {
      delete _ht_t[i];
    }
  }
  if (_isnew_ht_alpha_2) {
    delete _ht_alpha_2;
  }
  if (_isnew_ht_alpha_t) {
    delete _ht_alpha_t;
  }
  if (_isnew_ht_scale) {
    delete _ht_scale;
  }
  if(_spline_f2) {
    delete _spline_f2;
  }
  if(_spline_ft) {
    delete _spline_ft;
  }
}

void DIS_HT::update()
{
  std::vector<double> ht_x(_ht_x.size());
  std::vector<double> ht_f2(_ht_2.size());
  std::vector<double> ht_ft(_ht_t.size());
  for (size_t i = 0; i < ht_x.size(); i++)
  {
    ht_x[i] = *_ht_x[i];
    ht_f2[i] = *_ht_2[i];
    ht_ft[i] = *_ht_t[i];
  }
  _spline_ft->set_points(ht_x, ht_ft);
  _spline_f2->set_points(ht_x, ht_f2);
}

void DIS_HT::apply(const double q2, const double x, double& f2, double& fl) 
{
  static constexpr double q02 = 1.;
  double ft = f2 - fl;
  const double f2_cor = std::pow(x, *_ht_alpha_2) * (*_spline_f2)(x) * q02 / q2;
  const double ft_cor = std::pow(x, *_ht_alpha_t) * (*_spline_ft)(x) * q02 / q2;
  //printf("HT q2,x = %f,%f f2,ft = %f,%f\n", q2, x, f2_cor/f2, ft_cor/ft);
  f2 += f2_cor * (*_ht_scale);
  ft += ft_cor * (*_ht_scale);
  fl = f2 - ft;
}
