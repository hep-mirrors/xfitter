
/*
   @file ReactionFFABM_DISCC.cc
   @date 2017-10-09
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-10-09
*/

#include "ReactionFFABM_DISCC.h"
#include "ABM.h"
#include "DIS_HT.h"
#include "DIS_TMC.h"
#include "DIS_NUKE.h"
#include "xfitter_cpp_base.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include "cuba.h"
#include "cubature.h"
#include <spline.h>
struct integration_params_cuba {
  int intvar;
  double val;
  unsigned dataSetID;
  const BaseDISCC::ReactionData* rd;
  ReactionFFABM_DISCC* reaction;
  const double* br0;
  const double* br1;
  double mnucl;
};

// the class factories
extern "C" ReactionFFABM_DISCC* create() {
  return new ReactionFFABM_DISCC();
}

extern "C" {
double numufcalflux_(const double& e); // NOMAD E(nu) flux
double f2qcd_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2);
double flqcd_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2);
double f3qcd_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2);
double f2nucharm_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2, const int& nq);
double ftnucharm_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2, const int& nq);
double f3nucharm_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2, const int& nq);
// TODO to be removed
void sf_abkm_wrap_(const double& x, const double& q2,
                   const double& f2abkm, const double& flabkm, const double& f3abkm,
                   const double& f2cabkm, const double& flcabkm, const double& f3cabkm,
                   const double& f2babkm, const double& flbabkm, const double& f3babkm,
                   const int& ncflag, const double& charge, const double& polar,
                   const double& sin2thw, const double& cos2thw, const double& MZ, const int& nt=1);
void abkm_set_input_(const int& kschemepdfin, const int& kordpdfin,
                     const double& rmass8in, const double& rmass10in, const int& msbarmin,
                     double& hqscale1in, const double& hqscale2in, const int& flord);
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
  _FLAG_FAST = td->hasParam("_FLAG_FAST") ? td->getParamI("_FLAG_FAST") : 1;

  // scales mu^2 = scalea1 * Q^2 + scaleb1 * 4*m_h^2 (default scalea1 = scaleb1 = 1.0)
  _hqscale1in = 1.0;
  _hqscale2in = 1.0;
  if (td->hasParam("scalea1"))
    _hqscale1in = *td->getParamD("scalea1");
  if(td->hasParam("scaleb1"))
    _hqscale2in = *td->getParamD("scaleb1");
  abm::set_hq_scales(_hqscale1in, _hqscale2in);

  // pole or MCbar running mass treatment (default pole)
  _msbarmin = false;
  if(td->hasParam("runm"))
  _msbarmin = *td->getParamD("runm");

  // O(alpha_S) F_L = O(alpha_S) F_2 + ordfl (default ordfl = 1)
  _ordfl = 1;
  if(td->hasParam("ordfl"))
  _ordfl = td->getParamI("ordfl");
  
  // control x range (certain PDF sets have limited x_min, x_max)
  if(td->hasParam("xbmin"))
    abm::set_xbmin(*td->getParamD("xbmin"));
  if(td->hasParam("xbmax"))
    abm::set_xbmax(*td->getParamD("xbmax"));

    abm::initgridconst();

  // Take the 3-flavour scheme as a default
  _kschemepdfin = 0;

  // heavy quark masses
  _mcPtr = td->getParamD("mch");
  _mbPtr = td->getParamD("mbt");
  
  // CKM matrix
  _ckm.resize(9);
  _ckm[0] = td->getParamD("Vud");
  _ckm[1] = td->getParamD("Vus");
  _ckm[2] = td->getParamD("Vub");
  _ckm[3] = td->getParamD("Vcd");
  _ckm[4] = td->getParamD("Vcs");
  _ckm[5] = td->getParamD("Vcb");
  _ckm[6] = td->getParamD("Vtd");
  _ckm[7] = td->getParamD("Vts");
  _ckm[8] = td->getParamD("Vtb");

  printf("---------------------------------------------\n");
  printf("INFO from ReactionFFABM_DISCC:\n");
  printf("FF ABM running mass def? T(rue), (F)alse: %c\n", _msbarmin ? 'T' : 'F');
  printf("O(alpha_S) F_L - O(alpha_S) F2 = %d\n", _ordfl);
  printf("factorisation scale for heavy quarks  is set to sqrt(%f * Q^2 + %f * 4m_q^2\n", _hqscale1in, _hqscale2in);
  printf("---------------------------------------------\n");

  unsigned termID = td->id;
  _order[termID] = OrderMap(td->getParamS("Order")) - 1;
  _orderHQ[termID] = (td->hasParam("OrderHQ")) ? OrderMap(td->getParamS("OrderHQ")) - 1 : -1;
  auto nBins = td->getNbins();
  BaseDISCC::ReactionData *rd = (BaseDISCC::ReactionData *)td->reactionData;
  if(rd->_integrated)
    nBins = rd->_integrated->getBinValuesQ2()->size();
  _f2abm[termID].resize(nBins);
  _flabm[termID].resize(nBins);
  _f3abm[termID].resize(nBins);

  _mzPtr = td->getParamD("Mz");
  _sin2thwPtr = td->getParamD("sin2thW");
  _cos2thw = 1 - *_sin2thwPtr;
}

//
void ReactionFFABM_DISCC::atIteration() {
  Super::atIteration ();

  abm::set_hq_masses(*_mcPtr, *_mbPtr);
  abm::update_ckm_matrix(_ckm);

  //if (_ht) {
  //  _ht->update();
  //}
  for (auto ht : _ht) {
    if (ht.second) {
      ht.second->update();
    }
  }

  // Flag for internal arrays
  for ( auto ds : _dsIDs)  {
    (_f2abm[ds])[0] = -100.;
    (_flabm[ds])[0] = -100.;
    (_f3abm[ds])[0] = -100.;
  }
}

int ReactionFFABM_DISCC::integrate_nomad(const int* ndim, const cubareal* inp, const int *ncomp, cubareal* val, void *params) {
  //(void)ndim; // Unused parameter
  const integration_params_cuba& pars = *(integration_params_cuba*)params;
  static const double mnucl = pars.mnucl;
  static constexpr double emin = 6.;
  //static constexpr double emin = 0.;
  //static constexpr double emin = 0.1;
  static constexpr double emax = 300.;
  static const double xmin = 1./(2.*mnucl*emax);
  static constexpr double xmax = 0.99;
  static constexpr double q2min = 1.;
  static const double smax = 2 * mnucl * emax + mnucl * mnucl;
  double q2max = smax * xmax;
  double x(-1.), q2(-1.), y(-1.), e(-1.), s(-1.), factor(-1.);
  if(pars.intvar == 1) { // E
    throw 42;
    if (*ndim == 2) {
      e = pars.val;
      s = 2 * mnucl * e + mnucl * mnucl;
      q2max = s * xmax;
      x = xmin + (xmax - xmin) * inp[0];
      q2 = q2min + (q2max - q2min) * inp[1];
      factor = (q2max - q2min) * (xmax - xmin);
    }
    else if (*ndim == 3) {
      x = xmin + (xmax - xmin) * inp[0];
      q2 = q2min + (q2max - q2min) * inp[1];
      e = emin + (emax - emin) * inp[2];
      s = 2 * mnucl * e + mnucl * mnucl;
      factor = numufcalflux_(e) * (emax - emin) * (q2max - q2min) * (xmax - xmin);
    }
    y = (q2 + x * x * mnucl * mnucl) / s / x;
  }
  else if(pars.intvar == 9) { // test
    //printf("SZ dupa\n");
    //double x = 0.015;
    //double q2 = 15.0;
    double x = 1e-3 + (1e-2 - 1e-3) * inp[0];
    double q2 = 10 + (100 - 10) * inp[1];
    double f2sum(0.), flsum(0.), xf3sum(0.);
    pars.reaction->calc_point(q2, x, pars.dataSetID, pars.rd, f2sum, flsum, xf3sum);
    //printf("f2sum = %f\n", f2sum);
    //throw 1;
    val[0] = xf3sum;
    val[0] *= (100 - 10) * (1e-2 - 1e-3);
    return 0;
  }
  else if(pars.intvar == 7) { // test
    double f2sum(0.), flsum(0.), xf3sum(0.);
    if(1==2){
      double x = 0.15;
      double q2 = 15.0;
      double y = 0.5;
      pars.reaction->calc_point(q2, x, pars.dataSetID, pars.rd, f2sum, flsum, xf3sum);
      double yplus = 1.0 + (1.0 - y) * (1.0 - y);
      if (1) {
        yplus -= 2.0 * std::pow(mnucl * x * y, 2.0)/q2;
      }
      double yminus = 1.0 - (1.0 - y) * (1.0 - y);
      auto charge = pars.rd->_isBeamNu ? -1 * pars.rd->_charge : pars.rd->_charge; // swap charge for nu beam, see InitTerm()
      val[0] = 0.5 * (1 + charge * pars.rd->_polarisation) * (yplus * f2sum - charge * yminus * xf3sum - y * y * flsum);
      //val[0] = 2*((1.-y-std::pow((mnucl*x*y),2.0)/q2)*f2sum+y*y/2.*(f2sum-flsum)+(y-y*y/2.)*xf3sum);
      double e = q2/(2*mnucl*x*y);
      printf("SZ = %f %f %f %f %f %f %f %f\n", x, y, q2, mnucl, f2sum*2, flsum*2, xf3sum*2/x, val[0]);
      double br = *pars.br0 / (1 + *pars.br1 / e);
      val[0] *= br;
      return 0;
      throw 1;
    }
    double e = 0.15910e+02;
    double s=2*mnucl*e;
    double xmin1=q2min/s;
    x = xmin1 + (xmax - xmin1) * inp[0];
    double ymin = q2min/s/x;
    y = ymin + (1 - ymin) * inp[1];
    q2 = s * x * y;
    factor = 1.0;
    factor *= (xmax - xmin1) * (1 - ymin);
    factor *= x;
    factor *= e;
    factor /= x;
    pars.reaction->calc_point(q2, x, pars.dataSetID, pars.rd, f2sum, flsum, xf3sum);
    double yplus = 1.0 + (1.0 - y) * (1.0 - y);
    if (1) {
      yplus -= 2.0 * std::pow(mnucl * x * y, 2.0)/q2;
    }
    double yminus = 1.0 - (1.0 - y) * (1.0 - y);
    auto charge = pars.rd->_isBeamNu ? -1 * pars.rd->_charge : pars.rd->_charge; // swap charge for nu beam, see InitTerm()
    val[0] = 0.5 * (1 + charge * pars.rd->_polarisation) * (yplus * f2sum - charge * yminus * xf3sum - y * y * flsum);
    val[0] *= factor;
    double br = *pars.br0 / (1 + *pars.br1 / e);
    val[0] *= br;
    return 0;
  }
  else if(pars.intvar == 8) { // test
    double e = 0.15910e+02;
    s = 2 * mnucl * e + mnucl * mnucl;
    double smax1 = 2 * mnucl * e + mnucl * mnucl;
    q2max = smax1 * xmax;
    double xmin1 = 1./(2.*mnucl*e);
    x = xmin1 + inp[0] * (xmax - xmin1);
    y = inp[1];
    q2 = s * x * y;
    factor = 1.0;
    factor *= (xmax - xmin1);
    factor *= x;
    factor *= e;
    double f2sum(0.), flsum(0.), xf3sum(0.);
    pars.reaction->calc_point(q2, x, pars.dataSetID, pars.rd, f2sum, flsum, xf3sum);
    double yplus = 1.0 + (1.0 - y) * (1.0 - y);
    if (1) {
      yplus -= 2.0 * std::pow(mnucl * x * y, 2.0)/q2;
    }
    double yminus = 1.0 - (1.0 - y) * (1.0 - y);
    auto charge = pars.rd->_isBeamNu ? -1 * pars.rd->_charge : pars.rd->_charge; // swap charge for nu beam, see InitTerm()
    val[0] = 0.5 * (1 + charge * pars.rd->_polarisation) * (yplus * f2sum - charge * yminus * xf3sum - y * y * flsum);
    val[0] *= factor;
    double br = *pars.br0 / (1 + *pars.br1 / e);
    val[0] *= br;
    factor /= x;
    return 0;
  }
  else if(pars.intvar == 10) { // test
    printf("SZ dupa\n");
    const double pi = 3.1415926535897932384626433832795029;
    double x = inp[0] * pi/2.;
    double y = inp[1] * pi/2.;
    val[0] = sin(x)*cos(y);
    val[0] *= pi/2. * pi/2.;
    return 0;
  }
  else if(pars.intvar == 11) { // E
    if (*ndim == 2) {
      e = pars.val;
      s=2*mnucl*e;
      double xmin1=q2min/s;
      x = xmin1 + (xmax - xmin1) * inp[0];
      double ymin = q2min/s/x;
      y = ymin + (1 - ymin) * inp[1];
      q2 = s * x * y;
      factor = 1.0;
      factor *= (xmax - xmin1) * (1 - ymin);
      //factor *= x;
    }
  }
  else if(pars.intvar == 2) { // xbj
    throw 42;
    if (*ndim == 2) {
      e = emin + (emax - emin) * inp[0];
      q2 = q2min + (q2max - q2min) * inp[1];
      s = 2 * mnucl * e + mnucl * mnucl;
      x = pars.val;
      y = (q2 + x * x * mnucl * mnucl) / s / x;
      factor = numufcalflux_(e) * (q2max - q2min) * (emax - emin);
    }
    else if (*ndim == 3) {
      q2max = smax * xmax;
      e = emin + (emax - emin) * inp[0];
      q2 = q2min + (q2max - q2min) * inp[1];
      x = xmin + (xmax - xmin) * inp[2];
      s = 2 * mnucl * e + mnucl * mnucl;
      y = (q2 + x * x * mnucl * mnucl) / s / x;
      factor = numufcalflux_(e) * (q2max - q2min) * (emax - emin) * (xmax - xmin);
    }
  }
  else if(pars.intvar == 12) { // xbj
    if (*ndim == 2) {
      x = pars.val;
      double alim = q2min/x/2./mnucl;
      double emin1 = std::max(std::min(alim,emax),emin);
      e = emin1 + (emax - emin1) * inp[0];
      double ymin1 = alim/e;
      y = ymin1 + (1 - ymin1) * inp[1];
      //s = 2 * mnucl * e + mnucl * mnucl;
      s = 2 * mnucl * e;
      q2 = s * x * y;
      //factor = numufcalflux_(e) * (q2max - q2min) * (emax - emin);
      factor = numufcalflux_(e);
      factor *= (emax - emin1) * (1 - ymin1);
      factor *= e;
      //e = 10.;
    }
  }
  else if(pars.intvar == 3) { // sqrts
    //throw 42;
    if (*ndim == 2) {
      double sqrtshat = pars.val;
      double xmin = 1./(sqrtshat*sqrtshat);
      x = xmin + (xmax - xmin) * inp[0];
      q2 = sqrtshat * sqrtshat / (1. / x - 1);
      if (q2 <= q2min) {
        val[0] = 0.;
        return 0;
      }
      e = emin + (emax - emin) * inp[1];
      s = 2 * mnucl * e + mnucl * mnucl;
      y = (q2 + x * x * mnucl * mnucl) / s / x;
      factor = numufcalflux_(e) * (xmax - xmin) * (emax - emin);
    }
  }
  else if(pars.intvar == 13) { // sqrts
    if (*ndim == 2) {
      double sqrtshat = pars.val;
      double shat = sqrtshat*sqrtshat;
      double emin1 = std::min(std::max(emin,(q2min+shat)/2./mnucl),emax);
      e = emin1 + (emax - emin1) * inp[0];
      double spmax = 2*mnucl*e;
      double xmin1 = (1./(1+shat/q2min));
      double xmax1 = (1-shat/spmax);
      //printf("SZ e: %f - %f, x: %f - %f\n", emin1, emax, xmin1, xmax1);
      x = xmin1 + (xmax1 - xmin1) * inp[1];
      q2 = shat * x / (1 - x);
      if (q2 < (q2min-1e-6)) {
        printf("q2 = %f\n", q2);
        throw 43.0;
        val[0] = 0.;
        return 0;
      }
      y = q2 / (2 * e * mnucl * x);
      factor = numufcalflux_(e);
      //factor = 1.;
      factor *= 1. / (2 * mnucl * (1 - x));
      //factor *= 2 * sqrtshat;
      factor *= (xmax1 - xmin1) * (emax - emin1);
      //factor /= (xmax1 - xmin1) * (emax - emin1);
      //e = 10.;
    }
  }
  if (y <= 0. || y >= 1.) {
    val[0] = 0.;
    return 0;
  }
  double f2sum(0.), flsum(0.), xf3sum(0.);
  pars.reaction->calc_point(q2, x, pars.dataSetID, pars.rd, f2sum, flsum, xf3sum);
  double yplus = 1.0 + (1.0 - y) * (1.0 - y);
  if (1) {
    yplus -= 2.0 * std::pow(mnucl * x * y, 2.0)/q2;
  }
  double yminus = 1.0 - (1.0 - y) * (1.0 - y);
  auto charge = pars.rd->_isBeamNu ? -1 * pars.rd->_charge : pars.rd->_charge; // swap charge for nu beam, see InitTerm()
  val[0] = 0.5 * (1 + charge * pars.rd->_polarisation) * (yplus * f2sum - charge * yminus * xf3sum - y * y * flsum);
  const double pi = 3.1415926535897932384626433832795029;
  double MW = *pars.rd->Mw;
  double GF = pars.reaction->_Gf;
  double kinfactor = (pow(MW, 4.) / pow((q2 + MW * MW), 2)) * GF * GF / (2 * pi * x) * pars.reaction->_convfac;
  //val[0] *= kinfactor;
  if(1==2){
    val[0] /= kinfactor;
    factor *= e;
  }
  val[0] *= factor;
  if ( pars.rd->_dataFlav == BaseDISCC::dataFlav::c ) {
    double br = *pars.br0 / (1 + *pars.br1 / e);
    val[0] *= br;
  }
  else if (1==1) {
    if (pars.rd->_dataFlav == BaseDISCC::dataFlav::incl || pars.rd->_dataFlav == BaseDISCC::dataFlav::l) {
      val[0] *= std::pow((MW * MW / (q2 + MW * MW)), 2.);
    }
  }
  if (val[0] != val[0] || 1./val[0] == 0.) {
    val[0] = 0.;
  }
  return 0;
}

// Place calculations in one function, to optimize calls.
void ReactionFFABM_DISCC::calcF2FL(int dataSetID) {
  if ( (_f2abm[dataSetID][0]< -99.) )
  {
    auto td = _tdDS[dataSetID];
    td->actualizeWrappers();
    abm::pdffillgrid();
    BaseDISCC::ReactionData *rd = (BaseDISCC::ReactionData *)td->reactionData;
    if (td->hasParam("nomad")) {
      int intvar = td->getParamI("nomad");
      auto& nomad_var  = *GetBinValues(td, "nomad_var");
      for (size_t i=0; i<nomad_var.size(); i++) {
        calc_integral(intvar, nomad_var[i], dataSetID, rd, _f2abm[dataSetID][i]);
      }
    }
    else {
      auto *q2p  = GetBinValues(td,"Q2"), *xp  = GetBinValues(td,"x");
      auto q2 = *q2p, x = *xp;
      for (size_t i=0; i<xp->size(); i++) {
        calc_point(q2[i], x[i], dataSetID, rd, _f2abm[dataSetID][i], _flabm[dataSetID][i], _f3abm[dataSetID][i]);
      }
    }
  }
}

// Calculates one data point as integral over Q2,x,E and returns values f2, fl, f3
void ReactionFFABM_DISCC::calc_integral(const int intvar, const double val, const int dataSetID, const BaseDISCC::ReactionData *rd, double& xsec_out)
{
  if (1 && _nuke[dataSetID]) {
    printf("SZ1\n");fflush(stdout);
    double f2(1.), fl(1.), f3(1.);
    double ret = _nuke[dataSetID]->apply(0.1, 10., f2, fl, f3);
    printf("SZ2 ret = %f\n", ret);fflush(stdout);
  }
  integration_params_cuba pars;
  pars.dataSetID = dataSetID;
  pars.rd = rd;
  pars.reaction = this;
  pars.intvar = *_tdDS[dataSetID]->getParamD("nomad");
  pars.val = val;
  static const double* br0 = _tdDS[dataSetID]->getParamD("br_cmu_0");
  static const double* br1 = _tdDS[dataSetID]->getParamD("br_cmu_1");
  pars.br0 = br0;
  pars.br1 = br1;
  static const double mpr = *_tdDS[dataSetID]->getParamD("mpr");
  static const double mnt = *_tdDS[dataSetID]->getParamD("mnt");
  double mnucl = (mpr+mnt)/2.;
  pars.mnucl = mnucl;
  const int NDIM = intvar > 0 ? 2 : 3;
  const int NCOMP = 1;
  void* USERDATA = &pars;
  const int NVEC = 1;
  //const double EPSREL = 1e-5;
  const double EPSREL = 1e-3;
  //const double EPSREL = 1e-2;
  const double EPSABS = 0;
  const int FLAGS = 0;
  const int MINEVAL = 0;
  const int MAXEVAL = NDIM == 2 ? 500000 : 1000000;
  const int KEY = 0;
  const char* STATEFILE = nullptr;
  void* SPIN = nullptr;
  int nregions = 0;
  int neval = 0;
  int fail = 0;
  cubareal cuba_integral[1], cuba_error[1], prob[1];
  Cuhre(NDIM, NCOMP, integrate_nomad, USERDATA, NVEC, EPSREL, EPSABS, FLAGS, MINEVAL, MAXEVAL, KEY, STATEFILE, SPIN, &nregions, &neval, &fail, cuba_integral, cuba_error, prob);
  printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
  for(int comp = 0; comp < NCOMP; ++comp )
    printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n", (double)cuba_integral[comp], (double)cuba_error[comp], (double)prob[comp]);
    xsec_out = cuba_integral[0];
  printf("SZ NOMAD nomad_var[%d] = %6.2f -> xsec = %.4e +- %.4e\n", pars.intvar, pars.val, xsec_out, cuba_error[0]);
}

// Calculates one data point at (Q2,x) and returns values f2, fl, f3
void ReactionFFABM_DISCC::calc_point(const double q2, const double x, const int dataSetID, const BaseDISCC::ReactionData *rd, double& f2out, double& flout, double& f3out)
{
  static constexpr int ncflag = 0;
  f2out = flout = f3out = 0.;
  if (q2 < 1.0) return;
  double f2(0), f2b(0), f2c(0), fl(0), flc(0), flb(0), f3(0), f3b(0), f3c(0);
  if (_FLAG_FAST) {
    if (rd->_dataFlav == BaseDISCC::dataFlav::incl || rd->_dataFlav == BaseDISCC::dataFlav::l) {
      f2 = calc_point_strfun(BaseDISCC::dataType::f2, BaseDISCC::dataFlav::l, q2, x, dataSetID, -1, rd->_charge);
      fl = calc_point_strfun(BaseDISCC::dataType::fl, BaseDISCC::dataFlav::l, q2, x, dataSetID, -1, rd->_charge);
      f3 = calc_point_strfun(BaseDISCC::dataType::f3, BaseDISCC::dataFlav::l, q2, x, dataSetID, -1, rd->_charge);
    }
    if (rd->_dataFlav == BaseDISCC::dataFlav::incl || rd->_dataFlav == BaseDISCC::dataFlav::c) {
      f2c = calc_point_strfun(BaseDISCC::dataType::f2, BaseDISCC::dataFlav::c, q2, x, dataSetID, -1, rd->_charge);
      flc = calc_point_strfun(BaseDISCC::dataType::fl, BaseDISCC::dataFlav::c, q2, x, dataSetID, -1, rd->_charge, &f2c);
      f3c = calc_point_strfun(BaseDISCC::dataType::f3, BaseDISCC::dataFlav::c, q2, x, dataSetID, -1, rd->_charge);
    }
    if (rd->_dataFlav == BaseDISCC::dataFlav::incl || rd->_dataFlav == BaseDISCC::dataFlav::b) {
      f2b = calc_point_strfun(BaseDISCC::dataType::f2, BaseDISCC::dataFlav::b, q2, x, dataSetID, -1, rd->_charge);
      flb = calc_point_strfun(BaseDISCC::dataType::fl, BaseDISCC::dataFlav::b, q2, x, dataSetID, -1, rd->_charge);
      f3b = calc_point_strfun(BaseDISCC::dataType::f3, BaseDISCC::dataFlav::b, q2, x, dataSetID, -1, rd->_charge);
    }
    //printf("%f %f %f\n", f2, fl, f3);
  }
  else {
    abkm_set_input_(_kschemepdfin, _order[dataSetID], *_mcPtr, *_mbPtr, _msbarmin, _hqscale1in, _hqscale2in, _ordfl);
    sf_abkm_wrap_(x, q2, f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b, ncflag, rd->_charge, rd->_polarisation, *_sin2thwPtr, _cos2thw, *_mzPtr);
  }
  double f3out_bar = 0.;
  if(_nuke[dataSetID] && _nuke[dataSetID]->need_f3bar()) {
    // need F3bar for nuclear corrections and antineutrino
    // we can calculate it now, because HT and TMC (calculated later) do not apply to F3
    int charge_bar = -1 * rd->_charge;
    double f2_bar(0), f2b_bar(0), f2c_bar(0), fl_bar(0), flc_bar(0), flb_bar(0), f3_bar(0), f3c_bar(0), f3b_bar(0);
    if (_FLAG_FAST) {
      if (rd->_dataFlav == BaseDISCC::dataFlav::incl || rd->_dataFlav == BaseDISCC::dataFlav::l) {
        f3_bar = calc_point_strfun(BaseDISCC::dataType::f3, BaseDISCC::dataFlav::l, q2, x, dataSetID, -1, charge_bar);
      }
      if (rd->_dataFlav == BaseDISCC::dataFlav::incl || rd->_dataFlav == BaseDISCC::dataFlav::c) {
        f3c_bar = calc_point_strfun(BaseDISCC::dataType::f3, BaseDISCC::dataFlav::c, q2, x, dataSetID, -1, charge_bar);
      }
      if (rd->_dataFlav == BaseDISCC::dataFlav::incl || rd->_dataFlav == BaseDISCC::dataFlav::b) {
        f3b_bar = calc_point_strfun(BaseDISCC::dataType::f3, BaseDISCC::dataFlav::b, q2, x, dataSetID, -1, charge_bar);
      }
    }
    else {
      sf_abkm_wrap_(x, q2, f2_bar, fl_bar, f3_bar, f2c_bar, flc_bar, f3c_bar, f2b_bar, flb_bar, f3b_bar, ncflag, charge_bar, rd->_polarisation, *_sin2thwPtr, _cos2thw, *_mzPtr);
    }
    f3out_bar = x * combine_flavours(rd, f3_bar, f3c_bar, f3b_bar);
  }
  if (_tmc[dataSetID]) {
    if(((f2 != 0 || fl != 0 || f3 != 0) && _tmc[dataSetID]->getFlagL()) || ((f2c != 0 || flc != 0 || f3c != 0) && _tmc[dataSetID]->getFlagC()) || ((f2b != 0 || flb != 0 || f3b != 0) && _tmc[dataSetID]->getFlagB()) ) {
      const bool flag_fl = true;
      //const bool flag_f3 = false;
      const bool flag_f3 = true; //[not implemented]
      _tmc[dataSetID]->apply(f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b, flag_fl, flag_f3, q2, x, ncflag, rd->_charge, rd->_polarisation, _cos2thw, *_mzPtr, this);
    }
  }
  // HT is applied only to F2 and FL light flavour part
  //if(_flag_ht[dataSetID]) {
  //  _ht->apply(q2, x, f2, fl);
  //}
  if(_ht[dataSetID]) {
    if (f2 != 0 || fl != 0 || f3 != 0) {
      //printf("%f %f %f\n", f2, fl, f3);
      _ht[dataSetID]->apply(q2, x, f2, fl, f3);
    }
  }      
  f2out = combine_flavours(rd, f2, f2c, f2b);
  flout = combine_flavours(rd, fl, flc, flb);
  f3out = x * combine_flavours(rd, f3, f3c, f3b);
  // apply nuclear corrections to the sum of light+c+b because corrections for charm and non-charm (kint=4,5) are not implemented
  if (_nuke[dataSetID]) {
    _nuke[dataSetID]->apply(q2, x, f2out, flout, f3out, &f3out_bar);
  }
}

double ReactionFFABM_DISCC::calc_point_strfun(const BaseDISCC::dataType ftype, const BaseDISCC::dataFlav flav, const double q2, const double x, const int dataSetID, const int order, const int charge, const double* f2c) {
  if (flav == BaseDISCC::dataFlav::b) {
    return 0;
  }
  int orderALL = (order >= 0) ? order : _order[dataSetID];
  int orderHQ = (order >= 0) ? order : _orderHQ[dataSetID];
  int orderFL = (order >= 0) ? order : _ordfl;
  abm::set_scheme_and_order(_kschemepdfin, orderALL, _msbarmin, orderFL, orderHQ);
  static constexpr int nt = 1; // proton
  static constexpr int ni = 24; // CC
  const int nb = charge > 0 ? 6 : 7;
  switch (flav) {
    case BaseDISCC::dataFlav::l:
      switch (ftype) {
        case BaseDISCC::dataType::f2:
          return f2qcd_(nb, nt, ni, x, q2) / 2.;
        case BaseDISCC::dataType::fl:
          return flqcd_(nb, nt, ni, x, q2) / 2.;
        case BaseDISCC::dataType::f3:
          return f3qcd_(nb, nt, ni, x, q2) / 2.;
      }
    case BaseDISCC::dataFlav::c:
      double f2c_calc = 0.;
      switch (ftype) {
        case BaseDISCC::dataType::f2:
          return f2nucharm_(nb, nt, ni, x, q2, 8) / 2.;
        case BaseDISCC::dataType::fl:
          if (f2c) {
            f2c_calc = *f2c;
          }
          else {
            f2c_calc = f2nucharm_(nb, nt, ni, x, q2, 8) / 2.;
          }
          return f2c_calc - ftnucharm_(nb, nt, ni, x, q2, 8) / 2.;
        case BaseDISCC::dataType::f3:
          return f3nucharm_(nb, nt, ni, x, q2, 8) / 2.;
      }
  }
  hf_errlog(28022501, "F: Unsupported structure function type or flavour");
  return 0; // avoid warning
}

double ReactionFFABM_DISCC::combine_flavours(const BaseDISCC::ReactionData* rd, const double f, const double fc, const double fb)
{
  switch (rd->_dataFlav)
  {
    case BaseDISCC::dataFlav::incl:
      return f + fc + fb;
    case BaseDISCC::dataFlav::c:
      return fc;
    case BaseDISCC::dataFlav::b:
      return fb;
    case BaseDISCC::dataFlav::l:
      return f;
    default:
      hf_errlog(28022501, "F: Unsupported flavour");
      return 0.; // avoid warning
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
