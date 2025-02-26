#include "DIS_HT.h"
#include "TermData.h"
#include <cmath>
#include <sstream>

void DIS_HT::init(TermData* td) {
  _flag_ht = true;
  // for arrays, expect comma-separated strings
  // each item should be either double or parameter name
  auto read_array = [td](const std::string& parname) {
    std::istringstream ss(td->getParamS(parname));
    std::string token;
    int counter = 0;
    //std::vector<std::unique_ptr<double> > result;
    std::vector<const double* > result;
    while(std::getline(ss, token, ','))
    {
      if (td->hasParam(token)) {
        //result.push_back(unique_ptr<double>(new double(*td->getParamD(token))));
        result.push_back(td->getParamD(token));
        //string msg = "I: using higher twist parameter " + parname + "[" + std::to_string(counter) + "] = \"" + token + "\" as xfitter parameter";
        //hf_errlog(24051700+counter, msg);
      }
      else {
        try {
          //result.push_back(unique_ptr<double>(new double(stod(token))));
          //TODO: remember these created parameters and delete them, see also createConstantParameter() in xfitter_pars.cc
          result.push_back(new double(stod(token)));
          //string msg = "I: using higher twist parameter " + parname + "[" + std::to_string(counter) + "] = \"" + token + "\" as constant value";
          //hf_errlog(24051701+counter, msg);
        }
        catch (const std::invalid_argument&) {
          string msg = "I: using higher twist parameter " + parname + "[" + std::to_string(counter) + "] = \"" + token + "\" is not a parameter name or double";
          hf_errlog(24051702, msg);
        }
      }
      counter++;
    }
    return result;
  };
  auto read_double = [td](const std::string& parname) {
    std::istringstream ss(td->getParamS(parname));
    std::string token;
    //std::unique_ptr<double> result;
    const double* result;
    while(std::getline(ss, token, ','))
    {
      if (td->hasParam(token)) {
        //string msg = "I: using higher twist parameter " + parname + " = \"" + token + "\" as xfitter parameter";
        //hf_errlog(24051700, msg);
        //result = unique_ptr<double>(new double(*td->getParamD(token)));
        result = td->getParamD(token);
      }
      else {
        try {
          //string msg = "I: using higher twist parameter " + parname + " = \"" + token + "\" as constant value";
          //hf_errlog(24051701, msg);
          //result = unique_ptr<double>(new double(stod(token)));
          //TODO: remember these created parameters and delete them, see also createConstantParameter() in xfitter_pars.cc
          result = new double(stod(token));
        }
        catch (const std::invalid_argument&) {
          string msg = "I: using higher twist parameter " + parname + " = \"" + token + "\" is not a parameter name or double";
          hf_errlog(24051702, msg);
        }
      }
    }
    return result;
  };
  _ht_x = read_array("ht_x");
  _ht_2 = read_array("ht_vals_f2");
  _ht_t = read_array("ht_vals_ft");
  _ht_alpha_2 = read_double("ht_val_alpha_f2");
  _ht_alpha_t = read_double("ht_val_alpha_ft");
}

void DIS_HT::update()
{
  if (_flag_ht) {
    std::vector<double> ht_x(_ht_x.size());
    std::vector<double> ht_f2(_ht_x.size());
    std::vector<double> ht_ft(_ht_x.size());
    for (size_t i = 0; i < ht_x.size(); i++)
    {
      ht_x[i] = *_ht_x[i];
      ht_f2[i] = *_ht_2[i];
      ht_ft[i] = *_ht_t[i];
    }
    _spline_ft.set_points(ht_x, ht_ft);
    _spline_f2.set_points(ht_x, ht_f2);
  }
}

void DIS_HT::apply(const double q2, const double x, double& f2, double& fl) 
{
  if (_flag_ht) 
  {
    static constexpr double q02 = 1.;
    double ft = f2 - fl;
    const double f2_cor = std::pow(x, *_ht_alpha_2) * _spline_f2(x) * q02 / q2;
    const double ft_cor = std::pow(x, *_ht_alpha_t) * _spline_ft(x) * q02 / q2;
    //printf("HT q2,x = %f,%f f2,ft = %f,%f\n", q2, x, f2_cor/f2, ft_cor/ft);
    f2 += f2_cor;
    ft += ft_cor;
    fl = f2 - ft;
  }
}
