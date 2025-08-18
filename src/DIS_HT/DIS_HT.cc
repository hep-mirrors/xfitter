#include "DIS_HT.h"
#include "spline.h"
#include "TermData.h"
#include <cmath>
#include <sstream>
#include <iostream>

DIS_HT::DIS_HT(TermData* td) {
  _spline_f2 = nullptr;
  _spline_ft = nullptr;
  _spline_f3 = nullptr;
  // for arrays, expect comma-separated strings
  // each item should be either double or parameter name
  auto read_array = [td](const std::string& parname, std::vector<bool>& isnew) {
    isnew.clear();
    std::vector<const double* > result;
    if (!td->hasParam(parname)) {
      return result;
    }
    std::istringstream ss(td->getParamS(parname));
    std::string token;
    int counter = 0;
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
  auto read_double = [td](const std::string& parname, bool& isnew, const double def = -999.) {
    const double* result = nullptr;
    if (!td->hasParam(parname)) {
      if (def != -999.) {
        result = new double(def);
        isnew = true;
      }
      else {
        isnew = false;
      }
      return result;
    }
    std::istringstream ss(td->getParamS(parname));
    std::string token;
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
  _ht_3 = read_array("ht_vals_f3", _isnew_ht_3);
  _ht_alpha_2 = read_double("ht_val_alpha_f2", _isnew_ht_alpha_2, 0.0);
  _ht_alpha_t = read_double("ht_val_alpha_ft", _isnew_ht_alpha_t, 0.0);
  _ht_alpha_3 = read_double("ht_val_alpha_f3", _isnew_ht_alpha_3, 0.0);
  _ht_scale = read_double("ht_val_scale", _isnew_ht_scale, 1.0);
  //printf("_ht_? size = %lu %lu %lu\n", _ht_2.size(), _ht_t.size(), _ht_3.size());
  if (_ht_2.size()) _spline_f2 = new tk::spline();
  if (_ht_t.size()) _spline_ft = new tk::spline();
  if (_ht_3.size()) _spline_f3 = new tk::spline();
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
  for(size_t i = 0; i < _isnew_ht_3.size(); i++) {
    if(_isnew_ht_3[i]) {
      delete _ht_3[i];
    }
  }
  if (_isnew_ht_alpha_2) {
    delete _ht_alpha_2;
  }
  if (_isnew_ht_alpha_t) {
    delete _ht_alpha_t;
  }
  if (_isnew_ht_alpha_3) {
    delete _ht_alpha_3;
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
  if(_spline_f3) {
    delete _spline_f3;
  }
}

void DIS_HT::update()
{
  std::vector<double> ht_x(_ht_x.size());
  std::vector<double> ht_f2(_ht_2.size());
  std::vector<double> ht_ft(_ht_t.size());
  std::vector<double> ht_f3(_ht_3.size());
  for (size_t i = 0; i < ht_x.size(); i++)
  {
    if(ht_x.size()) ht_x[i] = *_ht_x[i];
    if(ht_f2.size()) ht_f2[i] = *_ht_2[i];
    if(ht_ft.size()) ht_ft[i] = *_ht_t[i];
    if(ht_f3.size()) ht_f3[i] = *_ht_3[i];
  }
  if(ht_f2.size()) _spline_f2->set_points(ht_x, ht_f2);
  if(ht_ft.size()) _spline_ft->set_points(ht_x, ht_ft);
  if(ht_f3.size()) _spline_f3->set_points(ht_x, ht_f3);
}

void DIS_HT::apply(const double q2, const double x, double& f2, double& fl, double& f3) 
{
  static constexpr double q02 = 1.;
  double ft = f2 - fl;
  //printf("spline: %p %p %p\n", _spline_f2, _spline_ft, _spline_f3);
  const double f2_cor = _spline_f2 ? (std::pow(x, *_ht_alpha_2) * (*_spline_f2)(x) * q02 / q2) : 0.0;
  const double ft_cor = _spline_ft ? (std::pow(x, *_ht_alpha_t) * (*_spline_ft)(x) * q02 / q2) : 0.0;
  //const double f3_cor = _spline_f3 ? (std::pow(x, *_ht_alpha_3) * (*_spline_f3)(x) * q02 / q2) : 0.0;
  const double f3_cor = 0.;
  //printf("HT q2,x = %f,%f f2,ft,f3 = %f,%f,%f\n", q2, x, f2_cor/f2, ft_cor/ft, f3_cor/f3);
  f2 += f2_cor * (*_ht_scale);
  ft += ft_cor * (*_ht_scale);
  fl = f2 - ft;
  //f3 += f3_cor * (*_ht_scale);
  f3 += f3_cor * (*_ht_scale) / x;
}
