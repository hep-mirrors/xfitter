
/*
   @file ReactionFFABM_DISCC.cc
   @date 2017-10-09
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-10-09
*/

#include "ReactionFFABM_DISCC.h"

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
    double rmassp[50];
    double rcharge[150];
  };
  extern COMMON_masses masses_;
}


// Initialize at the start of the computation
int ReactionFFABM_DISCC::atStart(const string &s)
{
  int isout = Super::atStart(s);

  // scales mu^2 = scalea1 * Q^2 + scaleb1 * 4*m_h^2 (default scalea1 = scaleb1 = 1.0)
  double hqscale1in = 1.0;
  double hqscale2in = 1.0;
  if(checkParam("scalea1"))
    hqscale1in = GetParam("scalea1");
  if(checkParam("scaleb1"))
    hqscale2in = GetParam("scaleb1");

  // pole or MCbar running mass treatment (default pole)
  bool msbarmin = false;
  if(checkParam("runm"))
    msbarmin = GetParamI("runm");

  // O(alpha_S) F_L = O(alpha_S) F_2 + ordfl (default ordfl = 1)
  int ordfl = 1;
  if(checkParam("ordfl"))
    ordfl = GetParamI("ordfl");

  initgridconst_();

  // Take the 3-flavour scheme as a default
  int kschemepdfin = 0;

  // heavy quark masses
  double rmass8in = GetParam("mch");
  masses_.rmass[7] = rmass8in;
  masses_.rcharge[7] = 0.6666666;
  _mc = rmass8in;
  double rmass10in = GetParam("mbt");
  masses_.rmass[9] = rmass10in;
  masses_.rcharge[9] = 0.3333333;
  _mb = rmass10in;

  printf("---------------------------------------------\n");
  printf("INFO from ABKM_init:\n");
  printf("FF ABM running mass def? T(rue), (F)alse: %c\n", msbarmin ? 'T' : 'F');
  printf("O(alpha_S) F_L - O(alpha_S) F2 = %d\n", ordfl);
  printf("---------------------------------------------\n");
  printf("factorisation scale for heavy quarks  is set to sqrt(%f * Q^2 + %f * 4m_q^2\n", hqscale1in, hqscale2in);

  const string order = GetParamS("Order");
  // NLO or NNLO: kordpdfin=1 NLO, kordpdfin=2 NNLO
  // this flag will set kordhq,kordalps,kordf2,kordfl,kordfl to same order
  const int kordpdfin = OrderMap( order) - 1;

  abkm_set_input_(kschemepdfin, kordpdfin, rmass8in, rmass10in, msbarmin, hqscale1in, hqscale2in, 1);
  abkm_set_input_orderfl_(ordfl);

  return isout;
}

void ReactionFFABM_DISCC::setDatasetParameters( int dataSetID, map<string,string> pars, map<string,double> parsDataset) {
  Super::setDatasetParameters(dataSetID, pars, parsDataset);
  // Allocate internal arrays:
  _f2abm[dataSetID].resize(GetNpoint(dataSetID));
  _flabm[dataSetID].resize(GetNpoint(dataSetID));
  _f3abm[dataSetID].resize(GetNpoint(dataSetID));
}

//
void ReactionFFABM_DISCC::initAtIteration() {

  Super::initAtIteration ();

  _mc = GetParam("mch");
  masses_.rmass[7] = _mc;
  _mb = GetParam("mbt");
  masses_.rmass[9] = _mb;

  //_asmz = alphaS(_mz);
  _mz = GetParam("Mz");
  _sin2thw = GetParam("sin2thW");
  _cos2thw = 1.0 - _sin2thw;

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
  if ( (_f2abm[dataSetID][0]< -99.) ) { // compute

      // CC
      int ncflag = 0;

      double charge = GetCharge(dataSetID);
      double polarity = GetPolarisation(dataSetID);

      // Get x,Q2 arrays:
      auto *q2p  = GetBinValues(dataSetID,"Q2"), *xp  = GetBinValues(dataSetID,"x");
      auto q2 = *q2p, x = *xp;

      // Number of data points
      const size_t Np = GetNpoint(dataSetID);

      double f2(0), f2b(0), f2c(0), fl(0), flc(0), flb(0), f3(0), f3b(0), f3c(0);

      for (size_t i=0; i<Np; i++) {
          if (q2[i]>1.0) {

              sf_abkm_wrap_(x[i], q2[i],
                           f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b,
                           ncflag, charge, polarity, _sin2thw, _cos2thw, _mz);
            }


          switch ( GetDataFlav(dataSetID) )
            {
            case dataFlav::incl :
              _f2abm[dataSetID][i] = f2 + f2c + f2b;
              _flabm[dataSetID][i] = fl + flc + flb;
              _f3abm[dataSetID][i] = x[i] * (f3 + f3c + f3b);
              break;
            case dataFlav::c :
              _f2abm[dataSetID][i] = f2c;
              _flabm[dataSetID][i] = flc;
              _f3abm[dataSetID][i] = x[i] * f3c;
              break;
            }
        }
  }
}

void ReactionFFABM_DISCC::F2 BASE_PARS
{
  calcF2FL(dataSetID);
  val = _f2abm[dataSetID];
}

void ReactionFFABM_DISCC::FL BASE_PARS
{
  calcF2FL(dataSetID);
  val = _flabm[dataSetID];
}

void ReactionFFABM_DISCC::xF3 BASE_PARS
{
  calcF2FL(dataSetID);
  val = _f3abm[dataSetID];
}
