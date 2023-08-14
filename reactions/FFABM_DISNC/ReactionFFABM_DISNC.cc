
/*
   @file ReactionFFABM_DISNC.cc
   @date 2017-09-29
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-09-29
*/

#include "ReactionFFABM_DISNC.h"
#include "xfitter_pars.h"
#include "xfitter_cpp_base.h"
#include <gsl/gsl_integration.h>
#include <spline.h>

// the class factories
extern "C" ReactionFFABM_DISNC *create()
{
  return new ReactionFFABM_DISNC();
}

// wrappers from:
//  ABM/src/sf_abkm_wrap.f
//  ABM/src/initgridconst.f
//  ABM/src/grid.f
extern "C"
{
void sf_abkm_wrap_(const double &x, const double &q2,
                   const double &f2abkm, const double &flabkm, const double &f3abkm,
                   const double &f2cabkm, const double &flcabkm, const double &f3cabkm,
                   const double &f2babkm, const double &flbabkm, const double &f3babkm,
                   const int &ncflag, const double &charge, const double &polar,
                   const double &sin2thw, const double &cos2thw, const double &MZ);
void abkm_set_input_(const int &kschemepdfin, const int &kordpdfin,
                     const double &rmass8in, const double &rmass10in, const int &msbarmin,
                     double &hqscale1in, const double &hqscale2in, const int &flagthinterface);
//void abkm_update_hq_masses_(const double& rmass8in, const double& rmass10in);
void abkm_set_input_orderfl_(const int &flord);
void initgridconst_();
void pdffillgrid_();

struct COMMON_masses
{
  double rmass[150];
  double rmassp[50];
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
void ReactionFFABM_DISNC::atStart()
{
  // do not call parent atStart(): it initialises QCDNUM
  // Super::atStart();
}

void ReactionFFABM_DISNC::initTerm(TermData *td)
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
  masses_.rcharge[7] = 0.6666666;
  _mbPtr = td->getParamD("mbt");
  masses_.rmass[9] = *_mbPtr;
  masses_.rcharge[9] = 0.3333333;

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
  if(_integrated.find(termID) != _integrated.end())
    nBins = _integrated[termID]->getBinValuesQ2()->size();
  _f2abm[termID].resize(nBins);
  _flabm[termID].resize(nBins);
  _f3abm[termID].resize(nBins);

  _mzPtr = td->getParamD("Mz");
  _sin2thwPtr = td->getParamD("sin2thW");

  // target mass correction
  if (td->hasParam("tmc"))
    _flag_tmc = *td->getParamD("tmc");
  else
    _flag_tmc = 0;
}

//
void ReactionFFABM_DISNC::atIteration()
{

  Super::atIteration();

  masses_.rmass[7] = *_mcPtr;
  masses_.rmass[9] = *_mbPtr;

  // need any TermData pointer to actualise PDFs and alpha_s
  // for the pdffillgrid_ call: use 1st one, this works properly
  // only if all terms have same evolution, decomposition etc.
  auto td = _tdDS.begin()->second;
  td->actualizeWrappers();
  pdffillgrid_();

  // Flag for internal arrays
  for (auto ds : _dsIDs)
  {
    (_f2abm[ds])[0] = -100.;
    (_flabm[ds])[0] = -100.;
    (_f3abm[ds])[0] = -100.;
  }
}

// Place calculations in one function, to optimize calls.
void ReactionFFABM_DISNC::calcF2FL(unsigned dataSetID)
{
  if ((_f2abm[dataSetID][0] < -99.))
  { // compute
    // use ref to termData:
    auto td = _tdDS[dataSetID];
    // NC
    int ncflag = 1;

    double charge = GetCharge(dataSetID);
    double polarity = GetPolarisation(dataSetID);

    // Get x,Q2 arrays:
    auto *q2p = GetBinValues(td, "Q2"), *xp = GetBinValues(td, "x");
    auto q2 = *q2p, x = *xp;

    // Number of data points
    // SZ perhaps GetNpoint does not work for integrated cross sections (not tested)
    //const size_t Np = GetNpoint(dataSetID);
    const size_t Np = xp->size();

    double f2(0), f2b(0), f2c(0), fl(0), flc(0), flb(0), f3(0), f3b(0), f3c(0);
    double cos2thw = 1.0 - *_sin2thwPtr;

    double maxdiff = -1.;
    for (size_t i = 0; i < Np; i++)
    {
      if (q2[i] > 1.0)
      {

        sf_abkm_wrap_(x[i], q2[i],
                      f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b,
                      ncflag, charge, polarity, *_sin2thwPtr, cos2thw, *_mzPtr);
        if(_flag_tmc) {
          // target mass corrections
          double mn = 0.938272;
          double gam = sqrt(1+4*x[i]*x[i]*mn*mn/q2[i]/q2[i]);
          double xi = 2*x[i]/(1+gam);
          auto integrate = [](double xip, void* params) {
            const integration_params& integrationParams = *(integration_params*)params;
            double f2(0), f2b(0), f2c(0), fl(0), flc(0), flb(0), f3(0), f3b(0), f3c(0);
            sf_abkm_wrap_(xip, integrationParams.q2[integrationParams.i],
                        f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b,
                        integrationParams.ncflag, integrationParams.charge, integrationParams.polarity, *integrationParams._sin2thwPtr, integrationParams.cos2thw, *integrationParams._mzPtr);
            return f2/xip/xip;
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
          // numerical integration
          gsl_function F;
          F.function = integrate;
          F.params = &pars;
          size_t alloc_space = 1000;
          gsl_integration_workspace * w = gsl_integration_workspace_alloc(alloc_space);
          double epsabs = 0;
          double epsrel = 1e-2;
          int key_param = 6;
          double result, error;
          gsl_integration_qag (&F, xi, 1.0, epsabs, epsrel, alloc_space, key_param, w, &result, &error);
          gsl_integration_workspace_free (w);
          // Simpson 3/8 integration
          double a = xi;
          //double b = log10(1/xi);
          double b = a*5;
          if (b > 0.999) 
            b = 0.999;
          double sim38 = (b-a)/8.*(integrate(a, &pars)+3*integrate((2*a+b)/3., &pars)+3*integrate((a+2*b)/3., &pars)+integrate(b, &pars));
          //printf("%f %f %f %f\n", integrate(a, &pars), integrate((2*a+b)/3., &pars), integrate((a+2*b)/3., &pars), integrate(b, &pars));
          double I = result;
          double f20 = f2;
          f2 = x[i]*x[i]/xi/xi/gam/gam/gam*f2 + 6*x[i]*x[i]*x[i]*mn*mn/q2[i]/gam/gam/gam/gam*I;
          double ft = f2 - fl;
          ft = x[i]*x[i]/xi/xi/gam*ft + 2*x[i]*x[i]*x[i]*mn*mn/q2[i]/gam/gam*I;
          fl = f2 - ft;
          double fl0 = fl;
          printf("SZ [x,q2 = %f %f] result +- error = %f +- %f [%f] sim38 = %f [%f] [%f]\n", x[i], q2[i], result, error, error/result, sim38, sim38/result-1, f2/f20-1);
          if (fabs(f20/f2-1) >  maxdiff)
            maxdiff = fabs(f20/f2-1);
        }
        if (_flag_ht[td->id]) {
          double q02 = 1.;
          double ft = f2 - fl;
          tk::spline spline_2;
          spline_2.set_points(_ht_x, _ht_2);
          //printf("%f", ft);
          f2 += pow(x[i], _ht_alpha_2) * spline_2(x[i]) * q02 / q2[i];
          tk::spline spline_t;
          spline_t.set_points(_ht_x, _ht_t);
          ft += pow(x[i], _ht_alpha_t) * spline_t(x[i]) * q02 / q2[i];
          fl = f2 - ft;
          //printf("  %f\n", ft);
        }
      }

      switch (GetDataFlav(dataSetID))
      {
        case dataFlav::incl:
          _f2abm[dataSetID][i] = f2 + f2c + f2b;
          _flabm[dataSetID][i] = fl + flc + flb;
          _f3abm[dataSetID][i] = x[i] * (f3 + f3c + f3b);
          break;
        case dataFlav::c:
          _f2abm[dataSetID][i] = f2c;
          _flabm[dataSetID][i] = flc;
          _f3abm[dataSetID][i] = x[i] * f3c;
          break;
        case dataFlav::b:
          _f2abm[dataSetID][i] = f2b;
          _flabm[dataSetID][i] = flb;
          _f3abm[dataSetID][i] = x[i] * f3b;
          break;
      }
    }
    printf("maxdiff = %f\n", maxdiff);
  }
}

void ReactionFFABM_DISNC::F2 BASE_PARS
{
  calcF2FL(td->id);
  val = _f2abm[td->id];
}

void ReactionFFABM_DISNC::FL BASE_PARS
{
  calcF2FL(td->id);
  val = _flabm[td->id];
}

void ReactionFFABM_DISNC::xF3 BASE_PARS
{
  calcF2FL(td->id);
  val = _f3abm[td->id];
}
