
/*
   @file ReactionFFABM_DISCC.cc
   @date 2017-10-09
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-10-09
*/

#include "ReactionFFABM_DISCC.h"
#include "xfitter_cpp_base.h"

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

// Place calculations in one function, to optimize calls.
void ReactionFFABM_DISCC::calcF2FL(int dataSetID) {
  if ( (_f2abm[dataSetID][0]< -99.) )
  { // compute
    // use ref to termData:
    auto td = _tdDS[dataSetID];
    BaseDISCC::ReactionData *rd = (BaseDISCC::ReactionData *)td->reactionData;

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
