
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
#include <cstring>

// the class factories
extern "C" ReactionFFABM_DISCC* create() {
  return new ReactionFFABM_DISCC();
}

extern "C" {
  double numufcalflux_(const double& e); // NOMAD E(nu) flux
}

extern "C" {
  double sd2_(double* acc, double (*f)(double*), void (*r)(int, double*, double*, double*));
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
  unsigned termID = td->id;

  // scales mu^2 = scalea1 * Q^2 + scaleb1 * 4*m_h^2 (default scalea1 = scaleb1 = 1.0)
  _hqscale1in = 1.0;
  _hqscale2in = 1.0;
  if (td->hasParam("scalea1"))
    _hqscale1in = *td->getParamD("scalea1");
  if(td->hasParam("scaleb1"))
    _hqscale2in = *td->getParamD("scaleb1");
  abm::set_hq_scales(_hqscale1in, _hqscale2in);

  // pole or MCbar running mass treatment (default pole)
  _msbarmin[termID] = false;
  if(td->hasParam("runm"))
  _msbarmin[termID] = *td->getParamD("runm");
  // O(alpha_S) F_L = O(alpha_S) F_2 + ordfl (default ordfl = 1)
  _ordfl[termID] = td->hasParam("ordfl") ? td->getParamI("ordfl") : 1;
  _order[termID] = OrderMap(td->getParamS("Order")) - 1;
  _orderHQ[termID] = (td->hasParam("OrderHQ")) ? OrderMap(td->getParamS("OrderHQ")) - 1 : -1;
  
  // control x range (certain PDF sets have limited x_min, x_max)
  if(td->hasParam("xbmin"))
    abm::set_xbmin(*td->getParamD("xbmin"));
  if(td->hasParam("xbmax"))
    abm::set_xbmax(*td->getParamD("xbmax"));
  
  abm::initgridconst();

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
  printf("FF ABM running mass def? T(rue), (F)alse: %c\n", _msbarmin[termID] ? 'T' : 'F');
  printf("Order = %d\n", _order[termID]);
  printf("Order HQ = %d\n", _orderHQ[termID]);
  printf("O(alpha_S) F_L - O(alpha_S) F2 = %d\n", _ordfl[termID]);
  printf("factorisation scale for heavy quarks  is set to sqrt(%f * Q^2 + %f * 4m_q^2\n", _hqscale1in, _hqscale2in);
  printf("---------------------------------------------\n");

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

  //auto td = _tdDS.begin()->second;
  //td->actualizeWrappers();
  //abm::pdffillgrid();
  _need_pdffillgrid = true;
}

int ReactionFFABM_DISCC::integrate_nomad_cubareal(const int* ndim, const cubareal* inp, const int *ncomp, cubareal* val, void *params) {
  double double_inp[*ndim];
  for(int i = 0; i < *ndim; i++) {
    double_inp[i] = inp[i];
  }
  double double_val = 0.;
  auto ret = integrate_nomad(ndim, double_inp, ncomp, &double_val, params);
  val[0] = double_val;
  return ret;
}

long unsigned ReactionFFABM_DISCC::_ncalls_nomad;

int ReactionFFABM_DISCC::integrate_nomad(const int* ndim, const cubareal* inp, const int *ncomp, cubareal* val, void *params) {
  _ncalls_nomad++;
  //(void)ndim; // Unused parameter
  const nomad_integration_params& pars = *(nomad_integration_params*)params;
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
  if(1==2 && pars.rd->_dataFlav == BaseDISCC::dataFlav::l) {
    double x = 0.2;
    double q2 = 3.0;
    //const int order = -1;
    double f2sum = 0.;
    double flsum = 0.;
    double xf3sum = 0.;
    //printf("%d\n", pars.rd->_dataFlav == BaseDISCC::dataFlav::l);
    pars.reaction->calc_point(q2, x, pars.dataSetID, pars.rd, f2sum, flsum, xf3sum);
    printf("F2 = %f\n", f2sum);
    throw 2;
  }
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
  if (pars.nomad_scalesemilepbr) {
    double br = *pars.br0 / (1 + *pars.br1 / e);
    val[0] *= br;
  }
  if (pars.nomad_scaleq2mw2) {
    //printf("MW1 = %f\n", MW);
    if (pars.rd->_dataFlav == BaseDISCC::dataFlav::incl || pars.rd->_dataFlav == BaseDISCC::dataFlav::l) {
      //printf("MW2 = %f\n", MW);
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
    bool need_update = td->actualizeWrappers();
    if (_need_pdffillgrid || need_update) {
      printf("ReactionFFABM_DISCC pdffillgrid()\n");
      abm::pdffillgrid();
      _need_pdffillgrid = false;
    }
    BaseDISCC::ReactionData *rd = (BaseDISCC::ReactionData *)td->reactionData;
    if (td->hasParam("nomad")) {
      // special case of integrated cross sections: TODO refactoring is needed
      int intvar = td->getParamI("nomad");
      int nomad_scaleq2mw2 = td->hasParam("nomad_scaleq2mw2") ? td->getParamI("nomad_scaleq2mw2") : 0;
      int nomad_scalesemilepbr = td->hasParam("nomad_scalesemilepbr") ? td->getParamI("nomad_scalesemilepbr") : 0;
      int nomad_verbose = td->hasParam("nomad_verbose") ? td->getParamI("nomad_verbose") : 0;
      double nomad_epsrel = td->hasParam("nomad_epsrel") ? *td->getParamD("nomad_epsrel") : 0.03; // somehow this provides accuracy ~ 0.1% for cuba
      int nomad_threads = td->hasParam("nomad_threads") ? td->getParamI("nomad_threads") : 0;
      std::string nomad_integrator = td->hasParam("nomad_integrator") ? td->getParamS("nomad_integrator") : "cuba";
      auto& nomad_var  = *GetBinValues(td, "nomad_var");
      for (size_t i=0; i<nomad_var.size(); i++) {
        _ncalls_nomad = 0;
        if (nomad_integrator == "cuba") {
          setenv("CUBACORES", (std::to_string(nomad_threads)).c_str(), 1);
          calc_integral_cuba(intvar, nomad_var[i], dataSetID, rd, _f2abm[dataSetID][i], nomad_scaleq2mw2, nomad_scalesemilepbr, nomad_epsrel, nomad_verbose);
        }
        else if (nomad_integrator == "sd2") {
          //calc_integral_sd2(intvar, nomad_var[i], dataSetID, rd, _f2abm[dataSetID][i], nomad_scaleq2mw2, nomad_scalesemilepbr, nomad_epsrel, nomad_verbose);
          _sd2_nomad_var = intvar;
          _sd2_nomad_pars = nomad_integration_params();
          _sd2_nomad_pars.dataSetID = dataSetID;
          _sd2_nomad_pars_static = make_integration_params(intvar, nomad_var[i], dataSetID, rd, nomad_scaleq2mw2, nomad_scalesemilepbr);
          double acc = 0.;
          //double s0=1.;
          double s0=sd2_(&acc, integrand_sd2, integrand_sd2_region);
          acc = s0 * nomad_epsrel;
          double s1=sd2_(&acc, integrand_sd2, integrand_sd2_region);
          if(nomad_verbose)
          {
            printf("s0 = %e\n", s0);
            printf("s1 = %e\n", s1);
          }
          _f2abm[dataSetID][i] = s1;
        }
        else {
          char str[256];
          sprintf(str, "F: unknown integrator %s", nomad_integrator.c_str());
          hf_errlog_(26061101, str, strlen(str));
        }
        if(nomad_verbose)
        {
          printf("_ncalls_nomad = %ld\n", _ncalls_nomad);
        }
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

double ReactionFFABM_DISCC::integrand_sd2(double x[]) {
  //integrate_nomad(const int* ndim, const cubareal* inp, const int *ncomp, cubareal* val, void *params)
  //return ret;
  //return 0.001;
  const int ndim = 2;
  const int ncomp = 1;
  cubareal val = 0.0;
  integrate_nomad(&ndim, x, &ncomp, &val, &_sd2_nomad_pars_static);
  return val;
}

void ReactionFFABM_DISCC::integrand_sd2_region(int ll, double* xx, double* aa, double* bb)
{
  (void)ll;
  (void)xx;
  const double del = 1e-8;
  //const double del = 0.0;
  aa[0] = del;
  bb[0] = 1.0 - del;
  aa[1] = del;
  bb[1] = 1.0 - del;
}

nomad_integration_params ReactionFFABM_DISCC::make_integration_params(const int intvar, const double val, const int dataSetID, const BaseDISCC::ReactionData *rd, const int nomad_scaleq2mw2, const int nomad_scalesemilepbr)
{
  nomad_integration_params pars;
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
  pars.nomad_scaleq2mw2 = nomad_scaleq2mw2;
  pars.nomad_scalesemilepbr = nomad_scalesemilepbr;
  return pars;
}

nomad_integration_params ReactionFFABM_DISCC::_sd2_nomad_pars_static;

// Calculates one data point as integral over Q2,x,E and returns values f2, fl, f3
void ReactionFFABM_DISCC::calc_integral_cuba(const int intvar, const double val, const int dataSetID, const BaseDISCC::ReactionData *rd, double& xsec_out, const int nomad_scaleq2mw2, const int nomad_scalesemilepbr, const double nomad_epsrel, const int nomad_verbose)
{
  // load nuclear correction tables once, otherwise they will be loaded multiple time in parallel
  //if (1 && _nuke[dataSetID]) {
  //  double f2(1.), fl(1.), f3(1.);
  //  double ret = _nuke[dataSetID]->apply(0.1, 10., f2, fl, f3);
  //}
  nomad_integration_params pars = make_integration_params(intvar, val, dataSetID, rd, nomad_scaleq2mw2, nomad_scalesemilepbr);
  const int NDIM = intvar > 0 ? 2 : 3;
  const int NCOMP = 1;
  void* USERDATA = &pars;
  const int NVEC = 1;
  const double EPSREL = nomad_epsrel;
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
  Cuhre(NDIM, NCOMP, integrate_nomad_cubareal, USERDATA, NVEC, EPSREL, EPSABS, FLAGS, MINEVAL, MAXEVAL, KEY, STATEFILE, SPIN, &nregions, &neval, &fail, cuba_integral, cuba_error, prob);
  if(nomad_verbose) printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
  for(int comp = 0; comp < NCOMP; ++comp )
    if(nomad_verbose) printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n", (double)cuba_integral[comp], (double)cuba_error[comp], (double)prob[comp]);
  xsec_out = cuba_integral[0];
  if(nomad_verbose) printf("SZ NOMAD nomad_var[%d] = %6.2f -> xsec = %.4e +- %.4e\n", pars.intvar, pars.val, xsec_out, cuba_error[0]);
}

// Calculates one data point at (Q2,x) and returns values f2, fl, f3
void ReactionFFABM_DISCC::calc_point(const double q2, const double x, const int dataSetID, const BaseDISCC::ReactionData *rd, double& f2out, double& flout, double& f3out)
{
  f2out = flout = f3out = 0.;
  if (q2 < 1.0) return;
  double f2(0), f2b(0), f2c(0), fl(0), flc(0), flb(0), f3(0), f3b(0), f3c(0);
  if (rd->_dataFlav == BaseDISCC::dataFlav::incl || rd->_dataFlav == BaseDISCC::dataFlav::l) {
    f2 = abm::calc_point_strfun_CC(abm::SFtype::f2, abm::SFflav::l, q2, x, -1, _order[dataSetID], _orderHQ[dataSetID], _ordfl[dataSetID], _msbarmin[dataSetID], rd->_charge);
    fl = abm::calc_point_strfun_CC(abm::SFtype::fl, abm::SFflav::l, q2, x, -1, _order[dataSetID], _orderHQ[dataSetID], _ordfl[dataSetID], _msbarmin[dataSetID], rd->_charge);
    f3 = abm::calc_point_strfun_CC(abm::SFtype::f3, abm::SFflav::l, q2, x, -1, _order[dataSetID], _orderHQ[dataSetID], _ordfl[dataSetID], _msbarmin[dataSetID], rd->_charge);
  }
  if (rd->_dataFlav == BaseDISCC::dataFlav::incl || rd->_dataFlav == BaseDISCC::dataFlav::c) {
    f2c = abm::calc_point_strfun_CC(abm::SFtype::f2, abm::SFflav::c, q2, x, -1, _order[dataSetID], _orderHQ[dataSetID], _ordfl[dataSetID], _msbarmin[dataSetID], rd->_charge);
    flc = abm::calc_point_strfun_CC(abm::SFtype::fl, abm::SFflav::c, q2, x, -1, _order[dataSetID], _orderHQ[dataSetID], _ordfl[dataSetID], _msbarmin[dataSetID], rd->_charge, &f2c);
    f3c = abm::calc_point_strfun_CC(abm::SFtype::f3, abm::SFflav::c, q2, x, -1, _order[dataSetID], _orderHQ[dataSetID], _ordfl[dataSetID], _msbarmin[dataSetID], rd->_charge);
  }
  if (rd->_dataFlav == BaseDISCC::dataFlav::incl || rd->_dataFlav == BaseDISCC::dataFlav::b) {
    f2b = abm::calc_point_strfun_CC(abm::SFtype::f2, abm::SFflav::b, q2, x, -1, _order[dataSetID], _orderHQ[dataSetID], _ordfl[dataSetID], _msbarmin[dataSetID], rd->_charge);
    flb = abm::calc_point_strfun_CC(abm::SFtype::fl, abm::SFflav::b, q2, x, -1, _order[dataSetID], _orderHQ[dataSetID], _ordfl[dataSetID], _msbarmin[dataSetID], rd->_charge);
    f3b = abm::calc_point_strfun_CC(abm::SFtype::f3, abm::SFflav::b, q2, x, -1, _order[dataSetID], _orderHQ[dataSetID], _ordfl[dataSetID], _msbarmin[dataSetID], rd->_charge);
  }
  double f3out_bar = 0.;
  if(_nuke[dataSetID] && _nuke[dataSetID]->need_f3bar()) {
    // need F3bar for nuclear corrections and antineutrino
    // we can calculate it now, because HT and TMC (calculated later) do not apply to F3
    int charge_bar = -1 * rd->_charge;
    double f2_bar(0), f2b_bar(0), f2c_bar(0), fl_bar(0), flc_bar(0), flb_bar(0), f3_bar(0), f3c_bar(0), f3b_bar(0);
    if (rd->_dataFlav == BaseDISCC::dataFlav::incl || rd->_dataFlav == BaseDISCC::dataFlav::l) {
      f3_bar = abm::calc_point_strfun_CC(abm::SFtype::f3, abm::SFflav::l, q2, x, -1, _order[dataSetID], _orderHQ[dataSetID], _ordfl[dataSetID], _msbarmin[dataSetID], charge_bar);
    }
    if (rd->_dataFlav == BaseDISCC::dataFlav::incl || rd->_dataFlav == BaseDISCC::dataFlav::c) {
      f3c_bar = abm::calc_point_strfun_CC(abm::SFtype::f3, abm::SFflav::c, q2, x, -1, _order[dataSetID], _orderHQ[dataSetID], _ordfl[dataSetID], _msbarmin[dataSetID], charge_bar);
    }
    if (rd->_dataFlav == BaseDISCC::dataFlav::incl || rd->_dataFlav == BaseDISCC::dataFlav::b) {
      f3b_bar = abm::calc_point_strfun_CC(abm::SFtype::f3, abm::SFflav::b, q2, x, -1, _order[dataSetID], _orderHQ[dataSetID], _ordfl[dataSetID], _msbarmin[dataSetID], charge_bar);
    }
    f3out_bar = x * combine_flavours(rd, f3_bar, f3c_bar, f3b_bar);
  }
  if (_tmc[dataSetID]) {
    if(((f2 != 0 || fl != 0 || f3 != 0) && _tmc[dataSetID]->getFlagL()) || ((f2c != 0 || flc != 0 || f3c != 0) && _tmc[dataSetID]->getFlagC()) || ((f2b != 0 || flb != 0 || f3b != 0) && _tmc[dataSetID]->getFlagB()) ) {
      static constexpr abm::SFproc ncflag = abm::SFproc::cc;
      _tmc[dataSetID]->apply(f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b, q2, x, ncflag, _order[dataSetID], _orderHQ[dataSetID], _ordfl[dataSetID], _msbarmin[dataSetID], rd->_charge, rd->_polarisation, _cos2thw, *_mzPtr);
    }
  }
  if(_ht[dataSetID]) {
    if (f2 != 0 || fl != 0 || f3 != 0) {
      // HT is applied only to F2 and FL light flavour part
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
