#include "DIS_TMC.h"
//#include "ReactionFFABM_DISNC.h"
//#include "ReactionFFABM_DISCC.h"
//#include "../reactions/FFABM_DISNC/ReactionFFABM_DISNC.h"
//#include "../reactions/FFABM_DISCC/ReactionFFABM_DISCC.h"
#include "TermData.h"
#include "hf_errlog.h"
#include <cmath>
#include <gsl/gsl_integration.h>

extern "C" {
  void sf_abkm_wrap_(const double& x, const double& q2,
                     const double& f2abkm, const double& flabkm, const double& f3abkm,
                     const double& f2cabkm, const double& flcabkm, const double& f3cabkm,
                     const double& f2babkm, const double& flbabkm, const double& f3babkm,
                     const int& ncflag, const double& charge, const double& polar,
                     const double& sin2thw, const double& cos2thw, const double& MZ, const int& nt=1);
  void sf_abkm_wrap_order_(const double &x, const double &q2,
                     const double &f2abkm, const double &flabkm, const double &f3abkm,
                     const double &f2cabkm, const double &flcabkm, const double &f3cabkm,
                     const double &f2babkm, const double &flbabkm, const double &f3babkm,
                     const int &ncflag, const double &charge, const double &polar,
                     const double &sin2thw, const double &cos2thw, const double &MZ, const int& kordpdfin, const int& nt=1);
}

struct integration_params {
  abm::SFproc sfproc;
  abm::SFtype sftype;
  abm::SFflav sfflav;
  double q2;
  int ncflag;
  int charge;
  double polarity;
  double cos2thw;
  double mz;
  int order;
  int orderDefault;
  int orderHQ;
  int orderFL;
  bool msbarm;
  //double xi;
  //int nt;
  // OBS
  int flag_calc_sf; // 0: F2, 1: FL, 2: F3
  int flag_calc_flav; // 0: l, 1: c, 2: b
  double calc_point_strfun(double xip) const {
    return abm::calc_point_strfun(sfproc, sftype, sfflav, q2, xip, order, orderDefault, orderHQ, orderFL, msbarm, charge, 1.-cos2thw, polarity, mz);
  }
};

DIS_TMC::DIS_TMC(TermData* td){
  const std::string& tmc_str = td->getParamS("tmc");
  _flag_l = _flag_c = _flag_b = false;
  if (tmc_str.find('l') != std::string::npos) {
    _flag_l = true;
  }
  if (tmc_str.find('c') != std::string::npos) {
    _flag_c = true;
  }
  if (tmc_str.find('b') != std::string::npos) {
    _flag_b = true;
  }
  _flag_f2 = td->hasParam("tmc_f2") ? td->getParamI("tmc_f2") : true;
  _flag_fl = td->hasParam("tmc_fl") ? td->getParamI("tmc_fl") : true;
  _flag_f3 = td->hasParam("tmc_f3") ? td->getParamI("tmc_f3") : true;
  if (td->hasParam("tmc_integration_method"))
  {
    auto s = td->getParamS("tmc_integration_method");
    if (s == "gsl")
      _integration_method = 0;
    else if (s == "simpson13")
      _integration_method = 1;
    else if (s == "power_simpson13")
      _integration_method = -1;
    else if (s == "simpson38")
      _integration_method = 2;
    else if (s == "boole")
      _integration_method = 3;
    //else if (s == "cuba")
    //  _integration_method = 4;
    else {
      string msg = "F: unknown tmc_integration_method = " + td->getParamS("tmc_integration_method");
      hf_errlog(24052201, msg);
    }
  }
  else
    _integration_method = 0; // "gsl"
  if (td->hasParam("tmc_xmin"))
    _xmin = *td->getParamD("tmc_xmin");
  else
    _xmin = 0.;
  if (td->hasParam("tmc_logxlogq2min"))
    _logxlogq2min = *td->getParamD("tmc_logxlogq2min");
  else
    _logxlogq2min = 0.;
  _mpr = td->getParamD("mpr");
}

void DIS_TMC::apply(double& f2, double& fl, double& f3, 
      double& f2c, double& flc, double& f3c, double& f2b, double& flb, double& f3b, 
      const double q2, const double x, const abm::SFproc ncflag, 
      const int orderDefault, const int orderHQ, const int orderFL, const bool msbarm,
      const int charge, const double polarity, const double cos2thw, const double mz) {
  if(!_flag_f2 && !_flag_fl && !_flag_f3) return;
  if ((_xmin == 0. || _xmin < x) && (_logxlogq2min == 0. || _logxlogq2min < log(x)*log(q2))) {
    if(_flag_l) {
      apply_one_flavour(f2, fl, f3, abm::SFflav::l, q2, x, ncflag, orderDefault, orderHQ, orderFL, msbarm, charge, polarity, cos2thw, mz);
    }
    if(_flag_c) {
      apply_one_flavour(f2c, flc, f3c, abm::SFflav::c, q2, x, ncflag, orderDefault, orderHQ, orderFL, msbarm, charge, polarity, cos2thw, mz);
    }
    if(_flag_b) {
      apply_one_flavour(f2b, flb, f3b, abm::SFflav::b, q2, x, ncflag, orderDefault, orderHQ, orderFL, msbarm, charge, polarity, cos2thw, mz);
    }
  }
}

double DIS_TMC::apply_one_flavour(double& f2, double& fl, double& f3, 
    const abm::SFflav flav, const double q2, const double x, const abm::SFproc ncflag, 
    const int orderDefault, const int orderHQ, const int orderFL, const bool msbarm, 
    const int charge, const double polarity, const double cos2thw, const double mz) {
  double mn = *_mpr;
  double gam = sqrt(1+4*x*x*mn*mn/q2);
  double xi = 2*x/(1+gam);
  if (xi>1) return 0.;
  auto integrate = [](double xip, void* params) {
    if(xip >= 1.) {
      return 0.;
    }
    const integration_params& integrationParams = *(integration_params*)params;
    return integrationParams.calc_point_strfun(xip)/xip/xip;
  };
  integration_params pars;
  pars.q2 = q2;
  pars.sfproc = ncflag;
  pars.sfflav = flav;
  pars.charge = charge;
  pars.polarity = polarity;
  pars.cos2thw = cos2thw;
  pars.mz = mz;
  pars.msbarm = msbarm;
  pars.order = 0; // compute integrals for TMC using LO structure functions
  pars.orderDefault = orderDefault;
  pars.orderHQ = orderHQ;
  pars.orderFL = orderFL;
  double I_F2 = 0.;
  double I_F3 = 0.;
  // gsl integration
  if (_integration_method == 0) {
    gsl_function F;
    F.function = integrate;
    F.params = &pars;
    size_t alloc_space = 1000;
    gsl_integration_workspace* w = nullptr;
    double epsabs = 0;
    double epsrel = 1e-3;
    int key_param = 6;
    double result, error;
    if(_flag_f2 || _flag_fl) {
      pars.sftype = abm::SFtype::f2;
      w = gsl_integration_workspace_alloc(alloc_space);
      gsl_integration_qag (&F, xi, 1.0, epsabs, epsrel, alloc_space, key_param, w, &result, &error);
      gsl_integration_workspace_free (w);
      I_F2 = result;
    }
    if (_flag_f3) {
      pars.sftype = abm::SFtype::f3;
      w = gsl_integration_workspace_alloc(alloc_space);
      gsl_integration_qag (&F, xi, 1.0, epsabs, epsrel, alloc_space, key_param, w, &result, &error);
      gsl_integration_workspace_free (w);
      I_F3 = result;
    }
  }
  // Simpson 1/3 integration
  else if (_integration_method == 1) {
    double a = xi;
    double b = a*5;
    if (b > 0.999) 
      b = 0.999;
    if(_flag_f2 || _flag_fl) {
      pars.sftype = abm::SFtype::f2;
      I_F2 = (b-a)/6.*(integrate(a, &pars)+4*integrate((a+b)/2., &pars)+integrate(b, &pars));
    }
    if (_flag_f3) {
      pars.sftype = abm::SFtype::f3;
      I_F3 = (b-a)/6.*(integrate(a, &pars)+4*integrate((a+b)/2., &pars)+integrate(b, &pars));
    }
  }
  // Simpson 1/3 integration with power approximation
  else if (_integration_method == -1) {
    double a = xi;
    double b = a*5;
    if (b > 0.999) 
      b = 0.999;
    if(_flag_f2 || _flag_fl) {
      pars.sftype = abm::SFtype::f2;
      double f1 = integrate(a, &pars);
      double alpha = log(f1)/log(xi);
      double f1pr = f1 - pow(a, alpha);
      double f2 = integrate((a+b)/2., &pars);
      double f2pr = f2 - pow((a+b)/2., alpha);
      double f3 = integrate(b, &pars);
      double f3pr = f3 - pow(b, alpha);
      double sim13 = (b-a)/6.*(f1pr+4*f2pr+f3pr);
      double part_power = (1-pow(xi, alpha+1))/(alpha+1);
      I_F2 = sim13 + part_power;
    }
    if (_flag_f3) {
      pars.sftype = abm::SFtype::f3;
      double f1 = integrate(a, &pars);
      double alpha = log(f1)/log(xi);
      double f1pr = f1 - pow(a, alpha);
      double f2 = integrate((a+b)/2., &pars);
      double f2pr = f2 - pow((a+b)/2., alpha);
      double f3 = integrate(b, &pars);
      double f3pr = f3 - pow(b, alpha);
      double sim13 = (b-a)/6.*(f1pr+4*f2pr+f3pr);
      double part_power = (1-pow(xi, alpha+1))/(alpha+1);
      I_F3 = sim13 + part_power;
    }
  }
  // Simpson 3/8 integration
  else if (_integration_method == 2) {
    double a = xi;
    //double b = log10(1/xi);
    double b = a*5;
    if (b > 0.999) 
      b = 0.999;
    if(_flag_f2 || _flag_fl) {
      pars.sftype = abm::SFtype::f2;
      I_F2 = (b-a)/8.*(integrate(a, &pars)+3*integrate((2*a+b)/3., &pars)+3*integrate((a+2*b)/3., &pars)+integrate(b, &pars));
    }
    if (_flag_f3) {
      pars.sftype = abm::SFtype::f3;
      I_F3 = (b-a)/8.*(integrate(a, &pars)+3*integrate((2*a+b)/3., &pars)+3*integrate((a+2*b)/3., &pars)+integrate(b, &pars));
    }
  }
  // Boole integration
  else if (_integration_method == 3) {
    double a = xi;
    //double b = log10(1/xi);
    //double b = a*5;
    //double maxfrac = (xi>0.25) ? 0.9 : 0.5;
    //double maxfrac = 0.5 + 0.4 * xi;
    double maxfrac = 0.4 + 0.6 * xi;
    //double maxfrac = sqrt(xi);
    double b = xi+(1-xi)*maxfrac;
    //if (b > 0.999) 
    //  b = 0.999;
    double h = (b-a)/4.;
    if(_flag_f2 || _flag_fl) {
      pars.sftype = abm::SFtype::f2;
      I_F2 = 2*h/45.*(7*integrate(a, &pars)+32*integrate(a+h, &pars)+12*integrate(a+2*h, &pars)+32*integrate(a+3*h, &pars)+integrate(b, &pars));
    }
    if (_flag_f3) {
      pars.sftype = abm::SFtype::f3;
      I_F3 = 2*h/45.*(7*integrate(a, &pars)+32*integrate(a+h, &pars)+12*integrate(a+2*h, &pars)+32*integrate(a+3*h, &pars)+integrate(b, &pars));
    }
  }
  // cuba integration
  /*else if (_integration_method == 4) {
    const int NDIM = 2;
    const int NCOMP = 1;
    pars.xi = xi;
    void* USERDATA = &pars;
    const int NVEC = 1;
    const double EPSREL = 1e-3;
    //const double EPSABS = 1e-12;
    const double EPSABS = 0;
    const int FLAGS = 0;
    const int MINEVAL = 0;
    const int MAXEVAL = 10000;
    const int KEY = 0;
    const char* STATEFILE = nullptr;
    void* SPIN = nullptr;
    int nregions = 0;
    int neval = 0;
    int fail = 0;
    cubareal cuba_integral[1], cuba_error[1], prob[1];
    Cuhre(NDIM, NCOMP, ReactionFFABM_DISNC::Integrand_Cuhre, USERDATA, NVEC, EPSREL, EPSABS, FLAGS, MINEVAL, MAXEVAL, KEY, STATEFILE, SPIN, &nregions, &neval, &fail, cuba_integral, cuba_error, prob);
    //printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    //nregions, neval, fail);
    //for(int comp = 0; comp < NCOMP; ++comp )
    //  printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
    //    (double)cuba_integral[comp], (double)cuba_error[comp], (double)prob[comp]);
    I = cuba_integral[0];
  }*/
  double f20 = f2;
  double fl0 = fl;
  double f30 = f3;
  pars.order = -1; // switch back to default order
  if(_flag_f2 || _flag_fl) {
    pars.sftype = abm::SFtype::f2;
    double f2_at_xi = integrate(xi, &pars)*xi*xi;
    double f2_tmc = x*x/xi/xi/gam/gam/gam*f2_at_xi + 6*x*x*x*mn*mn/q2/gam/gam/gam/gam*I_F2;
    pars.sftype = abm::SFtype::fl;
    double fl_at_xi = integrate(xi, &pars)*xi*xi;
    double ft_at_xi = f2_at_xi - fl_at_xi;
    double ft_tmc = x*x/xi/xi/gam*ft_at_xi + 2*x*x*x*mn*mn/q2/gam/gam*I_F2;
    double fl_tmc = f2_tmc - ft_tmc;
    if(_flag_f2) {
      f2 = f2_tmc;
    }
    if(_flag_fl) {
      fl = fl_tmc;
    }
  }
  if (_flag_f3) {
    pars.sftype = abm::SFtype::f3;
    double f3_at_xi = integrate(xi, &pars)*xi*xi;
    f3 = x/xi/gam/gam*f3_at_xi + 2*mn*mn/q2*x*x/gam/gam/gam*I_F3;
  }
  //printf("TMC [x,q2 = %.2f %.1f] F2,FL,F3: %.0f,%.0f,%.0f%%\n", x, q2, 100*(f2/f20-1), 100*(fl/fl0-1), 100*(f3/f30-1));
  //printf("SZ [x,q2 = %f %f] result +- error = %f +- %f [%f] sim38 = %f [%f] [%f]\n", x, q2, result, error, error/result, sim38, sim38/result-1, f2/f20-1);
  return f2/f20-1;
}
