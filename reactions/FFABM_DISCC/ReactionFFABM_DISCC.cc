
/*
   @file ReactionFFABM_DISCC.cc
   @date 2017-10-09
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-10-09
*/

#include "ReactionFFABM_DISCC.h"
#include "xfitter_cpp_base.h"

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include "cuba.h"
#include "cubature.h"
#include <spline.h>
struct integration_params_cuba {
  int intvar;
  double var;
  int ncflag;
  int charge;
  double polarity;
  double cos2thw;
  const double* _sin2thwPtr;
  const double* _mzPtr;
  BaseDISCC::dataFlav flav;
  ReactionFFABM_DISCC* reaction;
  unsigned dataSetID;
  BaseDISCC::ReactionData* rd;
  const double* br0;
  const double* br1;
  double mn;
};

// the class factories
extern "C" ReactionFFABM_DISCC* create() {
  return new ReactionFFABM_DISCC();
}

// TODO: in the future this reaction should be improved,
// currently it duplicates ReactionFFABM_DISNC with very few lines changed,
// common code should be separated

// wrappers from:
//  ABM/src/sf_abkm_wrap.f
//  ABM/src/initgridconst.f
//  ABM/src/grid.f
extern "C" {
void sf_abkm_wrap_(const double& x, const double& q2,
                   const double& f2abkm, const double& flabkm, const double& f3abkm,
                   const double& f2cabkm, const double& flcabkm, const double& f3cabkm,
                   const double& f2babkm, const double& flbabkm, const double& f3babkm,
                   const int& ncflag, const double& charge, const double& polar,
                   const double& sin2thw, const double& cos2thw, const double& MZ);
void sf_abkm_wrap_order_(const double &x, const double &q2,
                   const double &f2abkm, const double &flabkm, const double &f3abkm,
                   const double &f2cabkm, const double &flcabkm, const double &f3cabkm,
                   const double &f2babkm, const double &flbabkm, const double &f3babkm,
                   const int &ncflag, const double &charge, const double &polar,
                   const double &sin2thw, const double &cos2thw, const double &MZ, const int& kordpdfin);
double numufcalflux_(const double& e);
double nuke_fast_(const double& xb, const double& q2, const int& nsf, const int& ityp, const int& kint1, const int& kord, const int& ftyp, const float& syst);
void abkm_set_input_(const int& kschemepdfin, const int& kordpdfin,
                     const double& rmass8in, const double& rmass10in, const int& msbarmin,
                     double& hqscale1in, const double& hqscale2in, const int& flagthinterface);
//void abkm_update_hq_masses_(double& rmass8in, double& rmass10in);
void abkm_set_input_orderfl_(const int& flord);
void initgridconst_();
void pdffillgrid_();

struct COMMON_masses
{
  double rmass[150];
  double rmassp[150];
  double rcharge[150];
};
extern COMMON_masses masses_;

struct COMMON_gridset
{
  double delx1,delx2,delxp,dels1[8],dels2[8],xlog1,xlog2,x1,q2ini[8],q2min,q2max,xbmin,xbmax;
  int nxmgrid,nxpgrid,nspgrid,nsmgrid,khalf;
};
extern COMMON_gridset gridset_;

struct COMMON_constants_abkm
{
  double pi,alpha,alphady,rmpr,gfer2,sintc,sintw2,rmw,rmz,rgz,ckm[9],ckm2[9];
};
extern COMMON_constants_abkm constants_abkm_;
}


// Initialize at the start of the computation
void ReactionFFABM_DISCC::atStart()
{
  // do not call parent atStart(): it initialises QCDNUM
  // Super::atStart();
}

void ReactionFFABM_DISCC::initTerm(TermData *td)
{
  Super::initTerm(td);

  // scales mu^2 = scalea1 * Q^2 + scaleb1 * 4*m_h^2 (default scalea1 = scaleb1 = 1.0)
  double hqscale1in = 1.0;
  double hqscale2in = 1.0;
  if (td->hasParam("scalea1"))
    hqscale1in = *td->getParamD("scalea1");
  if(td->hasParam("scaleb1"))
    hqscale2in = *td->getParamD("scaleb1");

  // pole or MCbar running mass treatment (default pole)
  bool msbarmin = false;
  if(td->hasParam("runm"))
    msbarmin = *td->getParamD("runm");

  // O(alpha_S) F_L = O(alpha_S) F_2 + ordfl (default ordfl = 1)
  int ordfl = 1;
  if(td->hasParam("ordfl"))
    ordfl = td->getParamI("ordfl");
  
  // control x range (certain PDF sets have limited x_min, x_max)
  if(td->hasParam("xbmin"))
    gridset_.xbmin = *td->getParamD("xbmin");
  if(td->hasParam("xbmax"))
    gridset_.xbmax = *td->getParamD("xbmax");

  initgridconst_();

  // Take the 3-flavour scheme as a default
  int kschemepdfin = 0;

  // heavy quark masses
  _mcPtr = td->getParamD("mch");
  masses_.rmass[7] = *_mcPtr;
  _mbPtr = td->getParamD("mbt");
  masses_.rmass[9] = *_mbPtr;
  
  // CKM matrix
  constants_abkm_.ckm[0] = *td->getParamD("Vud");
  constants_abkm_.ckm[1] = *td->getParamD("Vus");
  constants_abkm_.ckm[2] = *td->getParamD("Vub");
  constants_abkm_.ckm[3] = *td->getParamD("Vcd");
  constants_abkm_.ckm[4] = *td->getParamD("Vcs");
  constants_abkm_.ckm[5] = *td->getParamD("Vcb");
  constants_abkm_.ckm[6] = *td->getParamD("Vtd");
  constants_abkm_.ckm[7] = *td->getParamD("Vts");
  constants_abkm_.ckm[8] = *td->getParamD("Vtb");

  printf("---------------------------------------------\n");
  printf("INFO from ABKM_init:\n");
  printf("FF ABM running mass def? T(rue), (F)alse: %c\n", msbarmin ? 'T' : 'F');
  printf("O(alpha_S) F_L - O(alpha_S) F2 = %d\n", ordfl);
  printf("---------------------------------------------\n");
  printf("factorisation scale for heavy quarks  is set to sqrt(%f * Q^2 + %f * 4m_q^2\n", hqscale1in, hqscale2in);

  const string order = td->getParamS("Order");
  // NLO or NNLO: kordpdfin=1 NLO, kordpdfin=2 NNLO
  // this flag will set kordhq,kordalps,kordf2,kordfl,kordfl to same order
  const int kordpdfin = OrderMap(order) - 1;

  abkm_set_input_(kschemepdfin, kordpdfin, *_mcPtr, *_mbPtr, msbarmin, hqscale1in, hqscale2in, 1);
  abkm_set_input_orderfl_(ordfl);

  unsigned termID = td->id;
  auto nBins = td->getNbins();
  BaseDISCC::ReactionData *rd = (BaseDISCC::ReactionData *)td->reactionData;
  if(rd->_integrated)
    nBins = rd->_integrated->getBinValuesQ2()->size();
  _f2abm[termID].resize(nBins);
  _flabm[termID].resize(nBins);
  _f3abm[termID].resize(nBins);

  _mzPtr = td->getParamD("Mz");
  _sin2thwPtr = td->getParamD("sin2thW");

  // target mass correction
  if (td->hasParam("tmc"))
  {
    _flag_tmc[termID] = td->getParamI("tmc");
    if (td->hasParam("tmc_integration_method"))
    {
      auto s = td->getParamS("tmc_integration_method");
      if (s == "gsl")
        _tmc_integration_method[termID] = 0;
      else if (s == "simpson13")
        _tmc_integration_method[termID] = 1;
      else if (s == "power_simpson13")
        _tmc_integration_method[termID] = -1;
      else if (s == "simpson38")
        _tmc_integration_method[termID] = 2;
      else if (s == "boole")
        _tmc_integration_method[termID] = 3;
      else if (s == "cuba")
        _tmc_integration_method[termID] = 4;
      else {
        string msg = "F: unknown tmc_integration_method = " + td->getParamS("tmc_integration_method");
        hf_errlog(24052201, msg);
      }
    }
    else
      _tmc_integration_method[termID] = 0; // "gsl"
    if (td->hasParam("tmc_xmin"))
      _tmc_xmin[termID] = *td->getParamD("tmc_xmin");
    else
      _tmc_xmin[termID] = 0.;
    if (td->hasParam("tmc_logxlogq2min"))
      _tmc_logxlogq2min[termID] = *td->getParamD("tmc_logxlogq2min");
    else
      _tmc_logxlogq2min[termID] = 0.;
    _tmc_mpr = td->getParamD("mpr");
  }
  else
    _flag_tmc[termID] = 0;
}

//
void ReactionFFABM_DISCC::atIteration() {

  Super::atIteration ();

  masses_.rmass[7] = *_mcPtr;
  masses_.rmass[9] = *_mbPtr;

  // need any TermData pointer to actualise PDFs and alpha_s
  // for the pdffillgrid_ call: use 1st one, this works properly
  // only if all terms have same evolution, decomposition etc.
  auto td = _tdDS.begin()->second;
  td->actualizeWrappers();
  pdffillgrid_();

  // Flag for internal arrays
  for ( auto ds : _dsIDs)  {
    (_f2abm[ds])[0] = -100.;
    (_flabm[ds])[0] = -100.;
    (_f3abm[ds])[0] = -100.;
  }

}

double ReactionFFABM_DISCC::apply_tmc(const int method, double& f2, double& fl, double& f3, const int flag_flavour, const std::valarray<double>& q2, const std::valarray<double>& x,
    const int ncflag, const int charge, const double polarity, const double cos2thw, const size_t i) {
  double mn = *_tmc_mpr;
  double gam = sqrt(1+4*x[i]*x[i]*mn*mn/q2[i]);
  double xi = 2*x[i]/(1+gam);
  if (xi>1) {throw 42;}
  auto integrate = [](double xip, void* params) {
    if(xip >= 1.) {
      return 0.;
    }
    const integration_params& integrationParams = *(integration_params*)params;
    double f2(0), f2b(0), f2c(0), fl(0), flc(0), flb(0), f3(0), f3b(0), f3c(0);
    if (integrationParams.order == -1)
      sf_abkm_wrap_(xip, integrationParams.q2[integrationParams.i],
                f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b,
                integrationParams.ncflag, integrationParams.charge, integrationParams.polarity, *integrationParams._sin2thwPtr, integrationParams.cos2thw, *integrationParams._mzPtr);
    else
      sf_abkm_wrap_order_(xip, integrationParams.q2[integrationParams.i],
                f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b,
                integrationParams.ncflag, integrationParams.charge, integrationParams.polarity, *integrationParams._sin2thwPtr, integrationParams.cos2thw, *integrationParams._mzPtr, integrationParams.order);
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
  };
  integration_params pars;
  pars.q2 = q2;
  pars.i = i;
  pars.ncflag = ncflag;
  pars.charge = charge;
  pars.polarity = polarity;
  pars.cos2thw = cos2thw;
  pars._sin2thwPtr = _sin2thwPtr;
  pars._mzPtr = _mzPtr;
  pars.flag_calc_fl = 0;
  pars.flag_flavour = flag_flavour;
  pars.order = 0;
  if (1==0) {
    double x0 = 0.01;
    double x1 = 1.;
    std::vector<double> q2s = {2., 5., 10., 100.};
    double q20 = 10.;
    printf("%10s%10s%15s%15s\n", "q2", "x", "f2", "fl");
    for (auto& q20 : q2s) {
      pars.q2 = q20;
      pars.flag_flavour = 1;
      for(int ii=0; ii<=100; ii++) {
        double xx = x0+ii*(x1-x0)/100.;
        pars.flag_calc_fl = 1;
        double fl = integrate(xx, &pars);
        pars.flag_calc_fl = 0;
        double f2 = integrate(xx, &pars);
        printf("%10.1f%10.6f%15.6e%15.6e\n", q20, xx, f2, fl);
      }
    }
    fflush(stdout);
    throw 42;
  }
  double I;
  // gsl integration
  if (method == 0) {
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
    //gsl_integration_qag (&F, x[i], 1.0, epsabs, epsrel, alloc_space, key_param, w, &result, &error);
    gsl_integration_workspace_free (w);
    I = result;
  }
  // Simpson 1/3 integration
  else if (method == 1) {
    double a = xi;
    double b = a*5;
    if (b > 0.999) 
      b = 0.999;
    double sim13 = (b-a)/6.*(integrate(a, &pars)+4*integrate((a+b)/2., &pars)+integrate(b, &pars));
    I = sim13;
  }
  // Simpson 1/3 integration with power approximation
  else if (method == -1) {
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
  else if (method == 2) {
    double a = xi;
    //double b = log10(1/xi);
    double b = a*5;
    if (b > 0.999) 
      b = 0.999;
    double sim38 = (b-a)/8.*(integrate(a, &pars)+3*integrate((2*a+b)/3., &pars)+3*integrate((a+2*b)/3., &pars)+integrate(b, &pars));
    I = sim38;
  }
  // Boole integration
  else if (method == 3) {
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
  else if (method == 4) {
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
  //f2 = x[i]*x[i]/xi/xi/gam/gam/gam*f2 + 6*x[i]*x[i]*x[i]*mn*mn/q2[i]/gam/gam/gam/gam*I;
  pars.order = -1;
  double f2_at_xi = integrate(xi, &pars)*xi*xi;
  //double f2_at_xi = integrate(xi, &pars)*xi;
  //printf("f2: %f %f\n", f2, f2_at_xi);
  double f2_tmc = x[i]*x[i]/xi/xi/gam/gam/gam*f2_at_xi + 6*x[i]*x[i]*x[i]*mn*mn/q2[i]/gam/gam/gam/gam*I;
  //double f2_orig = integrate(x[i], &pars)*x[i]*x[i];
  //f2 = f2 * f2_tmc / f2_orig;
  f2 = f2_tmc;
  //double ft = f2 - fl;
  //ft = x[i]*x[i]/xi/xi/gam*ft + 2*x[i]*x[i]*x[i]*mn*mn/q2[i]/gam/gam*I;
  pars.flag_calc_fl = 1;
  double fl_at_xi = integrate(xi, &pars)*xi*xi;
  //double fl_at_xi = integrate(xi, &pars)*xi;
  double ft_at_xi = f2_at_xi - fl_at_xi;
  double ft = x[i]*x[i]/xi/xi/gam*ft_at_xi + 2*x[i]*x[i]*x[i]*mn*mn/q2[i]/gam/gam*I;
  double fl_tmc = f2_tmc - ft;
  fl = fl_tmc;
  //fl = fl + x[i]*x[i]/gam/gam*(1-gam*gam)*f2_at_xi/xi/xi+mn*mn*x[i]*x[i]*x[i]/q2[i]/gam/gam/gam/gam*I;
  //double fl0 = fl;
  //pars.flag_calc_fl = 2;
  //double f3_at_xi = integrate(xi, &pars)*xi*xi;
  //f3 = x[i]/xi/gam/gam*f3_at_xi + 2*mn*mn/q2[i]*x[i]*x[i]/gam/gam/gam*I;
  /*double fl_orig = integrate(x[i], &pars)*x[i]*x[i];
  //fl = fl * fl_tmc / fl_orig;
  pars.order = -1;
  f2_at_xi = integrate(xi, &pars)*xi*xi;
  fl_at_xi = integrate(xi, &pars)*xi*xi;
  ft_at_xi = f2_at_xi - fl_at_xi;
  ft = x[i]*x[i]/xi/xi/gam*ft_at_xi + 2*x[i]*x[i]*x[i]*mn*mn/q2[i]/gam/gam*I;
  fl = f2-ft;*/
  //printf("SZ [x,q2 = %f %f] gsl = %f +- %f [%f] cuba = %f +- %f [%f] sim38 = %f [%f] TMC [%f]\n", x[i], q2[i], result, error, error/result, cuba_integral[0], cuba_error[0], cuba_integral[0]/result-1, sim38, sim38/result-1, f2/f20-1);
  //printf("SZ [x,q2 = %f %f] result +- error = %f +- %f [%f] sim38 = %f [%f] [%f]\n", x[i], q2[i], result, error, error/result, sim38, sim38/result-1, f2/f20-1);
  return f2/f20-1;
}

int ReactionFFABM_DISCC::Integrand_Cuhre(const int* ndim, const cubareal* inp, const int *ncomp, cubareal* val, void *params) {
  //(void)ndim; // Unused parameter
  const integration_params_cuba& integrationParams = *(integration_params_cuba*)params;
  ReactionFFABM_DISCC* reaction = integrationParams.reaction;
  unsigned dataSetID = integrationParams.dataSetID;
  //double mpr = 0.938272;
  double mpr = integrationParams.mn;
  //const double xmin = 0.;
  const double xmax = 1.0;
  //const double xmax = 0.80;
  const double xmin = 1e-4;
  //const double xmin = 1e-7;
  //const double xmax = 0.75;
  //const double q2min = 0.8;
  const double q2min = 1.;
  //const double q2max = 1000.;
  //const double xmin = 0.02;
  //const double xmax = 0.75;
  //const double q2min = 1.;
  //const double q2max = 20.;
  const double emin = 6.;
  const double emax = 300.;
  double x(-1.), q2(-1.), e(-1.), s(-1.), q2max(-1.);
  double flux(-1.);
  if(integrationParams.intvar == 1) { // E
    if (*ndim == 3) {
      double smax = 2*mpr*emax + mpr*mpr;
      q2max = smax*xmax;
      e = emin + (emax - emin) * inp[2];
      s = 2*mpr*e + mpr*mpr;
      flux = numufcalflux_(e);
      flux *= (emax - emin);
    }
    else if (*ndim == 2) {
      e = integrationParams.var;
      s = 2*mpr*e + mpr*mpr;
      q2max = s*xmax;
      flux = 1.;
    }
    x = xmin + (xmax - xmin) * inp[0];
    q2 = q2min + (q2max - q2min) * inp[1];
    flux *= (q2max - q2min) * (xmax - xmin);
  }
  else if(integrationParams.intvar == 2) { // xbj
    double extraf = 1.;
    if (*ndim == 3) {
      x = xmin + (xmax - xmin) * inp[2];
      extraf = (xmax - xmin);
    }
    else if (*ndim == 2) {
      x = integrationParams.var;
    }
    double smax = 2*mpr*emax + mpr*mpr;
    q2max = smax*xmax;
    e = emin + (emax - emin) * inp[0];
    q2 = q2min + (q2max - q2min) * inp[1];
    s = 2*mpr*e + mpr*mpr;
    flux = numufcalflux_(e);
    flux *= (q2max - q2min) * (emax - emin) * extraf;
  }
  else if(integrationParams.intvar == 3) { // sqrts
    double extraf = 1.;
    if (*ndim == 3) {
      double smax = 2*mpr*emax + mpr*mpr;
      q2max = smax*xmax;
      x = xmin + (xmax - xmin) * inp[0];
      q2 = q2min + (q2max - q2min) * inp[2];
      e = emin + (emax - emin) * inp[1];
      s = 2*mpr*e + mpr*mpr;
      extraf = (emax - emin);
    }
    else if (*ndim == 2) {
      //double sqrtshat = integrationParams.var;
      //e = emin + (emax - emin) * inp[0];
      //s = 2*mpr*e + mpr*mpr;
      //double smax = 2*mpr*emax + mpr*mpr;
      //q2max = smax*xmax;
      //q2 = q2min + (q2max - q2min) * inp[1];
      //x = sqrtshat / s;

      double sqrtshat = integrationParams.var;
      x = xmin + (xmax - xmin) * inp[0];
      e = (sqrtshat*sqrtshat-x*x*mpr*mpr)/(2*x*mpr);
      if (e < emin || e > emax) {
        val[0] = 0.;
        return 0;
      }
      double smax = 2*mpr*emax + mpr*mpr;
      q2max = smax*xmax;
      q2 = q2min + (q2max - q2min) * inp[1];
      s = 2*mpr*e + mpr*mpr;
      //s = sqrtshat / x;
      //e = (s - mpr*mpr)/ (mpr*2)
      //q2 = sqrtshat*sqrtshat*(1./x-1);
    }
    flux = numufcalflux_(e);
    flux *= (q2max - q2min) * (xmax - xmin) * extraf;
  }
  double y = (q2 + x*x*mpr*mpr) / s / x;
  if(integrationParams.intvar == 11) { // E dydx
    if (*ndim == 2) {
      //printf("SZ inp = %f,%f\n", inp[0], inp[1]);
      e = integrationParams.var;
      s = 2*mpr*e + mpr*mpr;
      q2max = s*xmax;
      flux = 1.;
      double ymin = 0.;
      double ymax = 1.;
      x = xmin + (xmax - xmin) * inp[0];
      y = ymin + (ymax - ymin) * inp[1];
      q2 = s * x * y - x * x * mpr * mpr;
      if (q2 <= q2min) {
        val[0] = 0.;
        return 0;
      }
      if (q2 > q2max) {
        throw 42;
      }
      flux *= (s-mpr*mpr) * x;
      flux *= (ymax - ymin) * (xmax - xmin);
    }
  }
  if (y <= 0. || y >= 1.) {
    val[0] = 0.;
    return 0;
  }
  double f2(0), f2b(0), f2c(0), fl(0), flc(0), flb(0), f3(0), f3b(0), f3c(0);
  sf_abkm_wrap_(x, q2,
                f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b,
                integrationParams.ncflag, integrationParams.charge, integrationParams.polarity, *integrationParams._sin2thwPtr, integrationParams.cos2thw, *integrationParams._mzPtr);
  //printf("SZ integral E,s = %f,%f q2,x,y = %f,%f,%f -> f2,fl,f3,f2c,flc,f3c = %f,%f,%f,%f,%f,%f integrationParams.ncflag = %d, integrationParams.charge = %d, integrationParams.polarity = %f, *integrationParams._sin2thwPtr = %f, integrationParams.cos2thw = %f, *integrationParams._mzPtr = %f\n", e, s, q2, x, y, f2,fl,f3,f2c,flc,f3c, integrationParams.ncflag, integrationParams.charge, integrationParams.polarity, *integrationParams._sin2thwPtr, integrationParams.cos2thw, *integrationParams._mzPtr);
  double f2sum(0),flsum(0),xf3sum(0);
  switch ( integrationParams.flav ) {
    case BaseDISCC::dataFlav::incl :
      if(reaction->_flag_tmc[dataSetID]) {
        if ((reaction->_tmc_xmin[dataSetID] == 0. || reaction->_tmc_xmin[dataSetID] < x) && (reaction->_tmc_logxlogq2min[dataSetID] == 0. || reaction->_tmc_logxlogq2min[dataSetID] < log(x)*log(q2))) {
          std::valarray<double> q2v = {q2};
          std::valarray<double> xv = {x};
          reaction->apply_tmc(reaction->_tmc_integration_method[dataSetID], f2, fl, f3, 1, q2v, xv, integrationParams.ncflag, integrationParams.charge, integrationParams.polarity, integrationParams.cos2thw, 0);
        }
      }
      if (reaction->_flag_ht[dataSetID]) {
        // for HT
        tk::spline spline_f2, spline_ft;
        std::vector<double> ht_x(reaction->_ht_x[dataSetID].size());
        std::vector<double> ht_f2(reaction->_ht_x[dataSetID].size());
        std::vector<double> ht_ft(reaction->_ht_x[dataSetID].size());
        if (reaction->_flag_ht[dataSetID]) {
          for (size_t i = 0; i < ht_x.size(); i++)
          {
            ht_x[i] = *reaction->_ht_x[dataSetID][i];
            ht_f2[i] = *reaction->_ht_2[dataSetID][i];
            ht_ft[i] = *reaction->_ht_t[dataSetID][i];
          }
          spline_ft.set_points(ht_x, ht_ft);
          spline_f2.set_points(ht_x, ht_f2);
        }
        //
        //printf("SZ HT1 f2,fl = %f,%f\n", f2,fl);
        double q02 = 1.;
        double ft = f2 - fl;
        f2 += std::pow(x, *reaction->_ht_alpha_2[dataSetID]) * spline_f2(x) * q02 / q2;
        tk::spline spline_t;
        spline_t.set_points(ht_x, ht_ft);
        ft += std::pow(x, *reaction->_ht_alpha_t[dataSetID]) * spline_ft(x) * q02 / q2;
        fl = f2 - ft;
        //printf("SZ HT2 f2,fl = %f,%f\n", f2,fl);
      }
      f2sum = f2 + f2c + f2b;
      flsum = fl + flc + flb;
      xf3sum = x * (f3 + f3c + f3b);
      break;
    case BaseDISCC::dataFlav::c :
      f2sum = f2c;
      flsum = flc;
      xf3sum = x * f3c;
      break;
  }
  double yplus = 1.0 + (1.0 - y) * (1.0 - y);
  double yminus = 1.0 - (1.0 - y) * (1.0 - y);
  if (integrationParams.charge > 0)
    val[0] = 0.5 * (1 + integrationParams.polarity) * (yplus * f2sum - yminus * xf3sum - y * y * flsum);
  else
    val[0] = 0.5 * (1 - integrationParams.polarity) * (yplus * f2sum + yminus * xf3sum - y * y * flsum);
  const double pi = 3.1415926535897932384626433832795029;
  double MW = *integrationParams.rd->Mw;
  double factor = (MW * MW * MW * MW / pow((q2 + MW * MW), 2)) * reaction->_Gf * reaction->_Gf / (2 * pi * x) * reaction->_convfac;
  //double factor = (MW * MW * MW * MW / pow((q2), 2)) * reaction->_Gf * reaction->_Gf / (2 * pi * x) * reaction->_convfac;
  val[0] *= factor;
  val[0] *= flux;
  if ( integrationParams.flav == BaseDISCC::dataFlav::c ) {
    double br = *integrationParams.br0 / (1 + *integrationParams.br1 / e);
    val[0] *= br;
  }
  if (val[0] != val[0] || 1./val[0] == 0.) val[0] = 0.;
  return 0;
}

void apply_nuke(const double x, const double q2, double& f2, double& fl, double& f3) {
  ;
}

// Place calculations in one function, to optimize calls.
void ReactionFFABM_DISCC::calcF2FL(int dataSetID) {
  if (0) {
    float syst = 0.;
    double xb=0.01;
    int nsf = 3;
    int ityp = 0;
    int kint1 = 2;
    int kord = 3;
    int ftyp = 4;
    for (int i = 0; i < 12; i++) {
      double q2 = 0.5 + 0.1 * i;
      double nuke = nuke_fast_(xb, q2, nsf, ityp, kint1, kord, ftyp, syst);
      printf("SZ xb,q2 = %f,%f -> nuke = %f\n", xb, q2, nuke);
    }
    throw 42;
  }
  if ( (_f2abm[dataSetID][0]< -99.) )
  //if ( 1 )
  { // compute
    // use ref to termData:
    auto td = _tdDS[dataSetID];
    BaseDISCC::ReactionData *rd = (BaseDISCC::ReactionData *)td->reactionData;

    //printf("SZ calcF2FL: dataSetID,rd->_dataFlav,ht = %d,%d,%d\n", dataSetID, rd->_dataFlav, _flag_ht[dataSetID]);

    // CC
    int ncflag = 0;

    // Get x,Q2 arrays:
    auto *q2p  = GetBinValues(td,"Q2"), *xp  = GetBinValues(td,"x");
    auto q2 = *q2p, x = *xp;

    // Number of data points
    // SZ getNbins does not work for integrated cross sections (returning number of bins)
    //const size_t Np = td->getNbins();
    const size_t Np = xp->size();

    double f2(0), f2b(0), f2c(0), fl(0), flc(0), flb(0), f3(0), f3b(0), f3c(0);
    double cos2thw = 1.0 - *_sin2thwPtr;

    for (size_t i=0; i<Np; i++) {
      if (q2[i]>1.0) {

        if (td->hasParam("nomad")) {
          //printf("SZ flux(E=%f) = %f\n", 20., numufcalflux_(20.));
          int nomad = Np;
          // PDFs
          td->actualizeWrappers();
          pdffillgrid_();
          integration_params_cuba pars;
          pars.ncflag = ncflag;
          pars.charge = rd->_charge;
          pars.polarity = rd->_polarisation;
          pars.cos2thw = cos2thw;
          pars._sin2thwPtr = _sin2thwPtr;
          pars._mzPtr = _mzPtr;
          pars.flav = rd->_dataFlav;
          pars.reaction = this;
          pars.dataSetID = dataSetID;
          pars.rd = rd;
          pars.intvar = *td->getParamD("nomad");
          pars.br0 =td->getParamD("br_cmu_0");
          pars.br1 =td->getParamD("br_cmu_1");
          auto& nomad_var  = *GetBinValues(td,"nomad_var");
          //double mpr = 0.938272;
          double mp = 0.93827208816;
          double mn = 0.93956542052;
          double mnucl = (mp+mn)/2.;
          //double mpr = 0.938272 * 56;
          //mnucl *= 56;
          pars.mn = mnucl;
          bool flag_total = td->hasParam("nomad_total");
          i = -10;
          for (size_t ii=0; ii<nomad; ii++) {
            pars.var = nomad_var[ii];
            double error(0.);
            const int NDIM = flag_total ? 3 : 2;
            const int NCOMP = 1;
            void* USERDATA = &pars;
            const int NVEC = 1;
            //const double EPSREL = 1e-3;
            const double EPSREL = 1e-2;
            //const double EPSABS = 1e-12;
            const double EPSABS = 0;
            const int FLAGS = 0;
            const int MINEVAL = 0;
            const int MAXEVAL = flag_total ? 500000 : 200000;
            const int KEY = 0;
            const char* STATEFILE = nullptr;
            void* SPIN = nullptr;
            int nregions = 0;
            int neval = 0;
            int fail = 0;
            cubareal cuba_integral[1], cuba_error[1], prob[1];
            Cuhre(NDIM, NCOMP, Integrand_Cuhre, USERDATA, NVEC, EPSREL, EPSABS, FLAGS, MINEVAL, MAXEVAL, KEY, STATEFILE, SPIN, &nregions, &neval, &fail, cuba_integral, cuba_error, prob);
            printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
            for(int comp = 0; comp < NCOMP; ++comp )
              printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n", (double)cuba_integral[comp], (double)cuba_error[comp], (double)prob[comp]);
            _f2abm[dataSetID][ii] = cuba_integral[0];
            error = cuba_error[0];
            printf("SZ NOMAD nomad_var[%d] = %6.2f -> xsec = %.4e +- %.4e\n", pars.intvar, pars.var, _f2abm[dataSetID][ii], error);
            if (flag_total) {
              for (size_t iii=0; iii<nomad; iii++) {
                _f2abm[dataSetID][iii] = cuba_integral[0];
              }
              return;
            }
            //gsl_function F;
            //double result = gsl_sf_bessel_J0(4.5);
          }
          return;
        }

        sf_abkm_wrap_(x[i], q2[i],
                      f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b,
                      ncflag, rd->_charge, rd->_polarisation, *_sin2thwPtr, cos2thw, *_mzPtr);
      }


      switch ( rd->_dataFlav )
      {
        case BaseDISCC::dataFlav::incl :
          _f2abm[dataSetID][i] = f2 + f2c + f2b;
          _flabm[dataSetID][i] = fl + flc + flb;
          _f3abm[dataSetID][i] = x[i] * (f3 + f3c + f3b);
          break;
        case BaseDISCC::dataFlav::c :
          _f2abm[dataSetID][i] = f2c;
          _flabm[dataSetID][i] = flc;
          _f3abm[dataSetID][i] = x[i] * f3c;
          break;
        case BaseDISCC::dataFlav::b:
          _f2abm[dataSetID][i] = 0.0;
          _flabm[dataSetID][i] = 0.0;
          _f3abm[dataSetID][i] = 0.0;
          break;
      }
      //apply_tmc(3, f2, fl, f3, 1, q2, x, ncflag, rd->_charge, rd->_polarisation, cos2thw, i);
      if(td->hasParam("nuke_kint") && td->getParamI("nuke_kint") != 0 && 1) {
        //const int ityp = td->getParamI("nuke_ityp");
        const int ityp = 0;
        const int kint = td->getParamI("nuke_kint");
        const int kord = OrderMap(td->getParamS("Order"));
        const int ftyp = td->getParamI("nuke_ftyp");
        const float syst = 0.;
        double cor_f1 = nuke_fast_(x[i], q2[i], 1, ityp, kint, kord, ftyp, syst);
        double cor_f2 = nuke_fast_(x[i], q2[i], 2, ityp, kint, kord, ftyp, syst);
        double cor_f3 = nuke_fast_(x[i], q2[i], 3, ityp, kint, kord, ftyp, syst);
        //double cor_f3 = nuke_fast_(x[i], q2[i], 3, ityp, abs(kint), kord, ftyp, syst);
        //if (kint<0) {
        //  double cor_f3_bar = nuke_fast_(x[i], q2[i], 3, ityp, -1*kint, kord, ftyp, syst);
        //  cor_f3 = cor_f3*2 - cor_f3_bar;
        //}
        //printf("SZnuke %d %d %d %d %f %f %f %f %f\n", ityp, kint, kord, ftyp, cor_f1, cor_f2, cor_f3, x[i], q2[i]);
        double f1 = (_f2abm[dataSetID][i] - _flabm[dataSetID][i]) / (2 * x[i]);
        f1 *= cor_f1;
        _f2abm[dataSetID][i] *= cor_f2;
        _flabm[dataSetID][i] = _f2abm[dataSetID][i] - 2 * x[i] * f1;
        _f3abm[dataSetID][i] *= cor_f3;
      }
    }
  }
}

valarray<double> ReactionFFABM_DISCC::F2(TermData *td)
{
  calcF2FL(td->id);
  return _f2abm[td->id];
}

valarray<double> ReactionFFABM_DISCC::FL(TermData *td)
{
  calcF2FL(td->id);
  return _flabm[td->id];
}

valarray<double> ReactionFFABM_DISCC::xF3(TermData *td)
{
  calcF2FL(td->id);
  return _f3abm[td->id];
}
