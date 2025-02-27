
/*
   @file ReactionFFABM_DISCC.cc
   @date 2017-10-09
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-10-09
*/

#include "ReactionFFABM_DISCC.h"
#include "DIS_HT.h"
#include "DIS_TMC.h"
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
                   const double& sin2thw, const double& cos2thw, const double& MZ, const int& nt=1);
void sf_abkm_wrap_order_(const double &x, const double &q2,
                   const double &f2abkm, const double &flabkm, const double &f3abkm,
                   const double &f2cabkm, const double &flcabkm, const double &f3cabkm,
                   const double &f2babkm, const double &flbabkm, const double &f3babkm,
                   const int &ncflag, const double &charge, const double &polar,
                   const double &sin2thw, const double &cos2thw, const double &MZ, const int& kordpdfin, const int& nt=1);
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
  _cos2thw = 1 - *_sin2thwPtr;
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

  if (_ht) {
    _ht->update();
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
  static const double emin = 6.;
  static const double emax = 300.;
  static const double xmin = 1./(2.*mnucl*emax);
  static const double xmax = 0.99;
  static const double q2min = 1.;
  double x(-1.), q2(-1.), y(-1.), e(-1.), s(-1.), q2max(-1.), factor(-1.);
  if(pars.intvar == 1) { // E
    if (*ndim == 2) {
      e = pars.val;
      s = 2 * mnucl * e + mnucl * mnucl;
      q2max = s * xmax;
      x = xmin + (xmax - xmin) * inp[0];
      q2 = q2min + (q2max - q2min) * inp[1];
      factor = (q2max - q2min) * (xmax - xmin);
    }
    else if (*ndim == 3) {
      double smax = 2 * mnucl * emax + mnucl * mnucl;
      q2max = smax * xmax;
      x = xmin + (xmax - xmin) * inp[0];
      q2 = q2min + (q2max - q2min) * inp[1];
      e = emin + (emax - emin) * inp[2];
      s = 2 * mnucl * e + mnucl * mnucl;
      factor = numufcalflux_(e) * (emax - emin) * (q2max - q2min) * (xmax - xmin);
    }
    y = (q2 + x * x * mnucl * mnucl) / s / x;
  }
  else if(pars.intvar == 2) { // xbj
    if (*ndim == 2) {
      double smax = 2 * mnucl * emax + mnucl * mnucl;
      q2max = smax * xmax;
      e = emin + (emax - emin) * inp[0];
      q2 = q2min + (q2max - q2min) * inp[1];
      s = 2 * mnucl * e + mnucl * mnucl;
      x = pars.val;
      y = (q2 + x * x * mnucl * mnucl) / s / x;
      factor = numufcalflux_(e) * (q2max - q2min) * (emax - emin);
    }
    else if (*ndim == 3) {
      double smax = 2 * mnucl * emax + mnucl * mnucl;
      q2max = smax * xmax;
      e = emin + (emax - emin) * inp[0];
      q2 = q2min + (q2max - q2min) * inp[1];
      x = xmin + (xmax - xmin) * inp[2];
      s = 2 * mnucl * e + mnucl * mnucl;
      y = (q2 + x * x * mnucl * mnucl) / s / x;
      factor = numufcalflux_(e) * (q2max - q2min) * (emax - emin) * (xmax - xmin);
    }
  }
  else if(pars.intvar == 3) { // sqrts
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
  if (y <= 0. || y >= 1.) {
    val[0] = 0.;
    return 0;
  }
  double f2sum(0.), flsum(0.), xf3sum(0.);
  pars.reaction->calc_point(q2, x, pars.dataSetID, pars.rd, f2sum, flsum, xf3sum);
  double yplus = 1.0 + (1.0 - y) * (1.0 - y);
  double yminus = 1.0 - (1.0 - y) * (1.0 - y);
  auto charge = pars.rd->_charge;
  if (pars.rd->_isBeamNu) {
    charge *= -1; // swap back, see InitTerm()
  }
  if (charge > 0)
    val[0] = 0.5 * (1 + pars.rd->_polarisation) * (yplus * f2sum - yminus * xf3sum - y * y * flsum);
  else
    val[0] = 0.5 * (1 - pars.rd->_polarisation) * (yplus * f2sum + yminus * xf3sum - y * y * flsum);
  const double pi = 3.1415926535897932384626433832795029;
  double MW = *pars.rd->Mw;
  double GF = pars.reaction->_Gf;
  double kinfactor = (pow(MW, 4.) / pow((q2 + MW * MW), 2)) * GF * GF / (2 * pi * x) * pars.reaction->_convfac;
  val[0] *= kinfactor;
  val[0] *= factor;
  if ( pars.rd->_dataFlav == BaseDISCC::dataFlav::c ) {
    double br = *pars.br0 / (1 + *pars.br1 / e);
    val[0] *= br;
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
    pdffillgrid_();
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
  if(1 && rd->_nuke_ftyp) {
    printf("SZ1\n");fflush(stdout);
    double ret = nuke_fast_(0.1, 10., 1, 0, 2, 3, 2, 0.);
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
  //const double EPSREL = 1e-3;
  const double EPSREL = 1e-2;
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
  static constexpr int nt = 1;
  f2out = 0.;
  flout = 0.;
  f3out = 0.;
  if (q2 < 1.0) return;
  double f2(0), f2b(0), f2c(0), fl(0), flc(0), flb(0), f3(0), f3b(0), f3c(0);
  sf_abkm_wrap_(x, q2, f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b, ncflag, rd->_charge, rd->_polarisation, *_sin2thwPtr, _cos2thw, *_mzPtr);
  bool calc_f3bar = rd->_nuke_ftyp && rd->_nuke_kint < 0; // for nuclear corrections
  //calc_f3bar = false;
  double f2_bar(0), f2b_bar(0), f2c_bar(0), fl_bar(0), flc_bar(0), flb_bar(0), f3_bar(0), f3c_bar(0), f3b_bar(0);
  if(calc_f3bar) {
    int charge_bar = -1 * rd->_charge;
    sf_abkm_wrap_(x, q2, f2_bar, fl_bar, f3_bar, f2c_bar, flc_bar, f3c_bar, f2b_bar, flb_bar, f3b_bar, ncflag, charge_bar, rd->_polarisation, *_sin2thwPtr, _cos2thw, *_mzPtr, nt);
  }
  if (_tmc[dataSetID]) {
    const bool flag_fl = true;
    const bool flag_f3 = false;
    //const bool flag_f3 = true; [not implemented]
    _tmc[dataSetID]->apply(f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b, flag_fl, flag_f3, q2, x, ncflag, rd->_charge, rd->_polarisation, _cos2thw, *_mzPtr);
  }
  if(_flag_ht[dataSetID]) {
    _ht->apply(q2, x, f2, fl);
  }
  combine_flavours(rd, f2, f2c, f2b, f2out);
  combine_flavours(rd, fl, flc, flb, flout);
  combine_flavours(rd, f3, f3c, f3b, f3out);
  f3out *= x;
  double f3out_bar(0.);
  if(calc_f3bar) {
    combine_flavours(rd, f3_bar, f3c_bar, f3b_bar, f3out_bar);
    f3out_bar *= x;
  }
  if(rd->_nuke_ftyp) 
  {
    static constexpr int ityp = 0;
    static constexpr float syst = 0.;
    double cor_f1 = nuke_fast_(x, q2, 1, ityp, rd->_nuke_kint, rd->_nuke_kord, rd->_nuke_ftyp, syst);
    double cor_f2 = nuke_fast_(x, q2, 2, ityp, rd->_nuke_kint, rd->_nuke_kord, rd->_nuke_ftyp, syst);
    double cor_f3 = nuke_fast_(x, q2, 3, ityp, rd->_nuke_kint, rd->_nuke_kord, rd->_nuke_ftyp, syst);
    if (rd->_nuke_kint < 0) 
    {
      int kint_bar = -1 * rd->_nuke_kint;
      double cor_f3_bar = nuke_fast_(x, q2, 3,  ityp, kint_bar, rd->_nuke_kord, rd->_nuke_ftyp, syst);
      cor_f3 = ((f3out + f3out_bar) * cor_f3 - f3out_bar * cor_f3_bar) / f3out;
    }
    double f1 = (f2out - flout) / (2 * x);
    f1 *= cor_f1;
    f2out *= cor_f2;
    flout = f2out - 2 * x * f1;
    f3out *= cor_f3;
  }
}

void ReactionFFABM_DISCC::combine_flavours(const BaseDISCC::ReactionData* rd, const double f, const double fc, const double fb, double& fout)
{
  switch (rd->_dataFlav)
  {
    case BaseDISCC::dataFlav::incl:
      fout = f + fc + fb;
      break;
    case BaseDISCC::dataFlav::c:
      fout = fc;
      break;
    case BaseDISCC::dataFlav::b:
      fout = fb;
      break;
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
