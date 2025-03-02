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
  double q2;
  int ncflag;
  int charge;
  double polarity;
  double cos2thw;
  double _mz;
  int flag_calc_fl;
  int flag_flavour;
  int order;
  double xi;
  int nt;
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
    else if (s == "cuba")
      _integration_method = 4;
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

void DIS_TMC::apply(double& f2, double& fl, double& f3, double& f2c, double& flc, double& f3c, double& f2b, double& flb, double& f3b, 
    const bool flag_fl, const bool flag_f3, const double q2, const double x, const int ncflag, const int charge, const double polarity, const double cos2thw, const double mz, void* ptr_reaction) {
  if ((_xmin == 0. || _xmin < x) && (_logxlogq2min == 0. || _logxlogq2min < log(x)*log(q2))) {
    if(_flag_l) {
      /*if (ncflag == 1) {
        ReactionFFABM_DISNC* reaction = (ReactionFFABM_DISNC*)ptr_reaction;
        f2 = reaction->calc_point_strfun(ReactionBaseDISNC::dataType::f2, ReactionBaseDISNC::dataFlav::l, q2[i], x[i], dataSetID, -1, GetCharge(dataSetID));
        fl = reaction->calc_point_strfun(ReactionBaseDISNC::dataType::fl, ReactionBaseDISNC::dataFlav::l, q2[i], x[i], dataSetID, -1, GetCharge(dataSetID));
        f3 = reaction->calc_point_strfun(ReactionBaseDISNC::dataType::f3, ReactionBaseDISNC::dataFlav::l, q2[i], x[i], dataSetID, -1, GetCharge(dataSetID));
      }*/
      apply_one_flavour(f2, fl, f3, flag_fl, flag_f3, 1, q2, x, ncflag, charge, polarity, cos2thw, mz);
    }
    if(_flag_c) {
      apply_one_flavour(f2c, flc, f3c, flag_fl, flag_f3, 2, q2, x, ncflag, charge, polarity, cos2thw, mz);
    }
    if(_flag_b) {
      apply_one_flavour(f2b, flb, f3b, flag_fl, flag_f3, 3, q2, x, ncflag, charge, polarity, cos2thw, mz);
    }
  }
}

double DIS_TMC::apply_one_flavour(double& f2, double& fl, double& f3, const bool flag_fl, const bool flag_f3, 
  const int flag_flavour, const double q2, const double x, const int ncflag, const int charge, const double polarity, const double cos2thw, const double mz) {
  //printf("apply_tmc q2,x = %f,%f\n", q2, x);
  double mn = *_mpr;
  double gam = sqrt(1+4*x*x*mn*mn/q2);
  double xi = 2*x/(1+gam);
  if (xi>1) {throw 42;}
  auto integrate = [&](double xip, void* params) {
    if(xip >= 1.) {
      return 0.;
    }
    const integration_params& integrationParams = *(integration_params*)params;
    //printf("integrate q2,xip = %f,%f\n", integrationParams.q2, xip);
    double f2(0), f2b(0), f2c(0), fl(0), flc(0), flb(0), f3(0), f3b(0), f3c(0);
    if (_FLAG_FAST) {
      if (integrationParams.flag_calc_fl == 0) {
        if (integrationParams.flag_flavour == 1)
          return calc_point_strfun(rd, BaseDISCC::dataType::f2, BaseDISCC::dataFlav::l, q2, x, dataSetID, -1, rd->_charge)/xip/xip;
        else if (integrationParams.flag_flavour == 2)
          return f2c/xip/xip;
        else if (integrationParams.flag_flavour == 3)
          return f2b/xip/xip;
        else
          return 0.;
      }
      else if (integrationParams.flag_calc_fl == 1) {
        if (integrationParams.flag_flavour == 1)
          return fl/xip/xip;
        else if (integrationParams.flag_flavour == 2)
          return flc/xip/xip;
        else if (integrationParams.flag_flavour == 3)
          return flb/xip/xip;
        else
          return 0.;
      }
      else if (integrationParams.flag_calc_fl == 2) {
        if (integrationParams.flag_flavour == 1)
          return f3/xip/xip;
        else if (integrationParams.flag_flavour == 2)
          return f3c/xip/xip;
        else if (integrationParams.flag_flavour == 3)
          return f3b/xip/xip;
        else
          return 0.;
      }
      throw 1;
    }
    else {
      if (integrationParams.order >= 0) {
        abkm_set_input_(_kschemepdfin, integrationParams.order, *_mcPtr, *_mbPtr, _msbarmin, _hqscale1in, _hqscale2in, _ordfl);
      }
      sf_abkm_wrap_(xip, integrationParams.q2, f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b, integrationParams.ncflag, integrationParams.charge, integrationParams.polarity, 1.-integrationParams.cos2thw, integrationParams.cos2thw, integrationParams._mz, integrationParams.order);
      if (integrationParams.order >= 0) {
        abkm_set_input_(_kschemepdfin, _order, *_mcPtr, *_mbPtr, _msbarmin, _hqscale1in, _hqscale2in, _ordfl);
      }
      if (integrationParams.flag_calc_fl == 0) {
        if (integrationParams.flag_flavour == 1)
          return f2/xip/xip;
        else if (integrationParams.flag_flavour == 2)
          return f2c/xip/xip;
        else if (integrationParams.flag_flavour == 3)
          return f2b/xip/xip;
        else
          return 0.;
      }
      else if (integrationParams.flag_calc_fl == 1) {
        if (integrationParams.flag_flavour == 1)
          return fl/xip/xip;
        else if (integrationParams.flag_flavour == 2)
          return flc/xip/xip;
        else if (integrationParams.flag_flavour == 3)
          return flb/xip/xip;
        else
          return 0.;
      }
      else if (integrationParams.flag_calc_fl == 2) {
        if (integrationParams.flag_flavour == 1)
          return f3/xip/xip;
        else if (integrationParams.flag_flavour == 2)
          return f3c/xip/xip;
        else if (integrationParams.flag_flavour == 3)
          return f3b/xip/xip;
        else
          return 0.;
      }
      throw 1;
    }
  };
  integration_params pars;
  pars.q2 = q2;
  pars.ncflag = ncflag;
  pars.charge = charge;
  pars.polarity = polarity;
  pars.cos2thw = cos2thw;
  pars._mz = mz;
  pars.flag_calc_fl = 0;
  pars.flag_flavour = flag_flavour;
  pars.order = 0;
  double I;
  // gsl integration
  if (_integration_method == 0) {
    gsl_function F;
    F.function = integrate;
    F.params = &pars;
    size_t alloc_space = 1000;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(alloc_space);
    double epsabs = 0;
    double epsrel = 1e-3;
    int key_param = 6;
    double result, error;
    gsl_integration_qag (&F, xi, 1.0, epsabs, epsrel, alloc_space, key_param, w, &result, &error);
    //gsl_integration_qag (&F, x, 1.0, epsabs, epsrel, alloc_space, key_param, w, &result, &error);
    gsl_integration_workspace_free (w);
    I = result;
  }
  // Simpson 1/3 integration
  else if (_integration_method == 1) {
    double a = xi;
    double b = a*5;
    if (b > 0.999) 
      b = 0.999;
    double sim13 = (b-a)/6.*(integrate(a, &pars)+4*integrate((a+b)/2., &pars)+integrate(b, &pars));
    I = sim13;
  }
  // Simpson 1/3 integration with power approximation
  else if (_integration_method == -1) {
    double a = xi;
    double b = a*5;
    if (b > 0.999) 
      b = 0.999;
    double f1 = integrate(a, &pars);
    double alpha = log(f1)/log(xi);
    double f1pr = f1 - pow(a, alpha);
    double f2 = integrate((a+b)/2., &pars);
    double f2pr = f2 - pow((a+b)/2., alpha);
    double f3 = integrate(b, &pars);
    double f3pr = f3 - pow(b, alpha);
    double sim13 = (b-a)/6.*(f1pr+4*f2pr+f3pr);
    double part_power = (1-pow(xi, alpha+1))/(alpha+1);
    I = sim13 + part_power;
  }
  // Simpson 3/8 integration
  else if (_integration_method == 2) {
    double a = xi;
    //double b = log10(1/xi);
    double b = a*5;
    if (b > 0.999) 
      b = 0.999;
    double sim38 = (b-a)/8.*(integrate(a, &pars)+3*integrate((2*a+b)/3., &pars)+3*integrate((a+2*b)/3., &pars)+integrate(b, &pars));
    I = sim38;
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
    double boole = 2*h/45.*(7*integrate(a, &pars)+32*integrate(a+h, &pars)+12*integrate(a+2*h, &pars)+32*integrate(a+3*h, &pars)+integrate(b, &pars));
    I = boole;
  }
  // cuba integration
  else if (_integration_method == 4) {
    /*const int NDIM = 2;
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
    I = cuba_integral[0];*/
  }
  //double I = result;
  //I = sim38;
  //I = 0;
  //I *= 5;
  double f20 = f2;
  //f2 = x*x/xi/xi/gam/gam/gam*f2 + 6*x*x*x*mn*mn/q2/gam/gam/gam/gam*I;
  pars.order = -1;
  double f2_at_xi = integrate(xi, &pars)*xi*xi;
  //double f2_at_xi = integrate(xi, &pars)*xi;
  //printf("f2: %f %f\n", f2, f2_at_xi);
  double f2_tmc = x*x/xi/xi/gam/gam/gam*f2_at_xi + 6*x*x*x*mn*mn/q2/gam/gam/gam/gam*I;
  //double f2_orig = integrate(x, &pars)*x*x;
  //f2 = f2 * f2_tmc / f2_orig;
  f2 = f2_tmc;
  if (flag_fl) {
    //double ft = f2 - fl;
    //ft = x*x/xi/xi/gam*ft + 2*x*x*x*mn*mn/q2/gam/gam*I;
    pars.flag_calc_fl = 1;
    double fl_at_xi = integrate(xi, &pars)*xi*xi;
    //double fl_at_xi = integrate(xi, &pars)*xi;
    double ft_at_xi = f2_at_xi - fl_at_xi;
    double ft = x*x/xi/xi/gam*ft_at_xi + 2*x*x*x*mn*mn/q2/gam/gam*I;
    double fl_tmc = f2_tmc - ft;
    fl = fl_tmc;
    //fl = fl + x*x/gam/gam*(1-gam*gam)*f2_at_xi/xi/xi+mn*mn*x*x*x/q2/gam/gam/gam/gam*I;
    //double fl0 = fl;
  }
  /*if (flag_f3) {
    pars.flag_calc_fl = 2;
    double f3_at_xi = integrate(xi, &pars)*xi*xi;
    f3 = x/xi/gam/gam*f3_at_xi + 2*mn*mn/q2*x*x/gam/gam/gam*I;
  }*/
  /*double fl_orig = integrate(x, &pars)*x*x;
  //fl = fl * fl_tmc / fl_orig;
  pars.order = -1;
  f2_at_xi = integrate(xi, &pars)*xi*xi;
  fl_at_xi = integrate(xi, &pars)*xi*xi;
  ft_at_xi = f2_at_xi - fl_at_xi;
  ft = x*x/xi/xi/gam*ft_at_xi + 2*x*x*x*mn*mn/q2/gam/gam*I;
  fl = f2-ft;*/
  //printf("SZ [x,q2 = %f %f] gsl = %f +- %f [%f] cuba = %f +- %f [%f] sim38 = %f [%f] TMC [%f]\n", x, q2, result, error, error/result, cuba_integral[0], cuba_error[0], cuba_integral[0]/result-1, sim38, sim38/result-1, f2/f20-1);
  //printf("SZ [x,q2 = %f %f] result +- error = %f +- %f [%f] sim38 = %f [%f] [%f]\n", x, q2, result, error, error/result, sim38, sim38/result-1, f2/f20-1);
  return f2/f20-1;
}
