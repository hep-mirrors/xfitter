
/*
   @file ReactionFFABM_DISNC.cc
   @date 2017-09-29
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-09-29
*/

#include "ReactionFFABM_DISNC.h"
#include "DIS_HT.h"
#include "DIS_TMC.h"
#include "DIS_NUKE.h"
#include "xfitter_pars.h"
#include "xfitter_cpp_base.h"
#include <gsl/gsl_integration.h>
#include <spline.h>
//#include "cuba.h"
//#include "cubature.h"
#include <unistd.h>
#include <sys/shm.h>
#include <sys/wait.h>
#include "xfitter_steer.h"


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
                   const double &sin2thw, const double &cos2thw, const double &MZ, const int &nt=1);
void sf_abkm_wrap_order_(const double &x, const double &q2,
                   const double &f2abkm, const double &flabkm, const double &f3abkm,
                   const double &f2cabkm, const double &flcabkm, const double &f3cabkm,
                   const double &f2babkm, const double &flbabkm, const double &f3babkm,
                   const int &ncflag, const double &charge, const double &polar,
                   const double &sin2thw, const double &cos2thw, const double &MZ, const int& kordpdfin, const int &nt=1);
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

  // parallel
  _ncpu[td->id] = 1;
  //printf("td->hasParam(threads) = %d\n", td->hasParam("threads"));
  if (td->hasParam("threads")) 
    _ncpu[td->id] = td->getParamI("threads");
  if (_ncpu[td->id] == -1) {
    _ncpu[td->id] = sysconf(_SC_NPROCESSORS_ONLN);
    hf_errlog(2023061401,"I: Will use "+std::to_string(_ncpu[td->id])+" threads");
  }
  //printf("_ncpu = %d\n", _ncpu[td->id]);
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

  if (_ht) {
    _ht->update();
  }

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
    auto td = _tdDS[dataSetID];
    td->actualizeWrappers();
    pdffillgrid_();

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

    // Calculate i-th data point, return values f2, fl, f3
    auto calc_point = [&](int i, double& f2out, double& flout, double& f3out) {
      f2out = flout = f3out = 0.;
      if (q2[i] > 1.0)
      {
        //printf("q2,x = %f,%f\n", q2[i], x[i]);
        sf_abkm_wrap_(x[i], q2[i],
                      f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b,
                      ncflag, charge, polarity, *_sin2thwPtr, cos2thw, *_mzPtr);
        double f3out_bar = 0.;
        if(_nuke[dataSetID] && _nuke[dataSetID]->need_f3bar()) {
          // need F3bar for nuclear corrections and antineutrino
          // we can calculate it now, because HT and TMC (calculated later) do not apply to F3
          int charge_bar = -1 * charge;
          double f2_bar(0), f2b_bar(0), f2c_bar(0), fl_bar(0), flc_bar(0), flb_bar(0), f3_bar(0), f3c_bar(0), f3b_bar(0);
          sf_abkm_wrap_(x[i], q2[i], f2_bar, fl_bar, f3_bar, f2c_bar, flc_bar, f3c_bar, f2b_bar, flb_bar, f3b_bar, ncflag, charge_bar, polarity, *_sin2thwPtr, cos2thw, *_mzPtr);
          f3out_bar = x[i] * combine_flavours(GetDataFlav(dataSetID), f3_bar, f3c_bar, f3b_bar);
        }
        if (_tmc[dataSetID]) {
          const bool flag_fl = true;
          const bool flag_f3 = false;
          //const bool flag_f3 = true; [not implemented]
          _tmc[dataSetID]->apply(f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b, flag_fl, flag_f3, q2[i], x[i], ncflag, charge, polarity, cos2thw, *_mzPtr);
        }
        if(_flag_ht[dataSetID]) {
          _ht->apply(q2[i], x[i], f2, fl);
        }      
        f2out = combine_flavours(GetDataFlav(dataSetID), f2, f2c, f2b);
        flout = combine_flavours(GetDataFlav(dataSetID), fl, flc, flb);
        f3out = x[i] * combine_flavours(GetDataFlav(dataSetID), f3, f3c, f3b);
            // apply nuclear corrections to the sum of light+c+b because corrections for charm and non-charm (kint=4,5) are not implemented
        if (_nuke[dataSetID]) {
          _nuke[dataSetID]->apply(q2[i], x[i], f2out, flout, f3out, &f3out_bar);
        }
      }
    };

    //printf("FFABM _ncpu = %d\n", _ncpu[dataSetID]);
    int ncpu =  xfitter::xf_ncpu(_ncpu[dataSetID]);
    //printf("FFABM ncpu = %d\n", ncpu);


    if (ncpu == 1) {
      for (size_t i = 0; i < Np; i++) {
        double f2, fl, f3;
        calc_point(i, f2, fl, f3);
        _f2abm[dataSetID][i] = f2;
        _flabm[dataSetID][i] = fl;
        _f3abm[dataSetID][i] = f3;
      }
    }
    else {
      // Shared memory for predictions
      int shmid;
      double* sharedArray;
      shmid = shmget(IPC_PRIVATE, sizeof(double) * Np * 3, IPC_CREAT | 0666);
      if (shmid < 0) {
        hf_errlog(2023060200,"F: Failed to create shared memory segment");
      }
      sharedArray = static_cast<double*>(shmat(shmid, nullptr, 0));
      if (sharedArray == reinterpret_cast<double*>(-1)) {
        hf_errlog(2023060201,"F: Failed to attach shared memory segment");
      }
      // define Chunks
      int chunkSize = Np / ncpu;
      int reminder  = Np % ncpu; 
      int first = 0;
      int startIndex = 0;
      int endIndex = 0;
      // loop over all
      for (int icpu = 0; icpu < std::min(ncpu, int(Np)); icpu++) {
        startIndex = endIndex;
        endIndex   = startIndex + chunkSize;
        if (icpu < reminder) {
          endIndex += 1;
        }
        pid_t pid = xfitter::xf_fork( std::min(ncpu, int(Np))  );
        if ( pid == 0) {       
          // close all open files (e.g. minuit.out.txt) to avoid multiple buffered output
          int fdlimit = (int)sysconf(_SC_OPEN_MAX);
          for (int i = STDERR_FILENO + 1; i < fdlimit; i++) {
            close(i);
          }
          //printf("CPU %d will compute %d -- %d\n", icpu, first+startIndex, first+endIndex-1);
          for (int i = first+startIndex; i < first+endIndex; i++) {
            //printf("CPU %d computing %d\n", icpu, i);
            double f2, fl, f3;
            calc_point(i, f2, fl, f3);
            //printf("CPU %d computing %d f2,fl,f3 = %f %f %f\n", icpu, i, f2,fl,f3);
            sharedArray[i] = f2;
            sharedArray[i+Np] = fl;
            sharedArray[i+2*Np] = f3;
          }
          exit(0);	    
        }
        else if (pid<0) {
          hf_errlog(2023060204,"F: Failed to create a fork process");	
        }
      }	
      // Wait ...
      int status;
      while (wait(&status) > 0);    
      //struct rusage usage;
      //while (wait3(&status, 0, &usage) > 0);    
      // Store result
      for (size_t i = 0; i<Np; i++) {
        _f2abm[dataSetID][i] = sharedArray[i];
        _flabm[dataSetID][i] = sharedArray[i+Np];
        _f3abm[dataSetID][i] = sharedArray[i+2*Np];
      }    
      // Detach and remove shared memory segments
      shmdt(sharedArray);
      shmctl(shmid, IPC_RMID, NULL);
    }
  }
}

double ReactionFFABM_DISNC::combine_flavours(const ReactionBaseDISNC::dataFlav flav, const double f, const double fc, const double fb)
{
  switch (flav)
  {
    case ReactionBaseDISNC::dataFlav::incl:
      return f + fc + fb;
    case ReactionBaseDISNC::dataFlav::c:
      return fc;
    case ReactionBaseDISNC::dataFlav::b:
      return fb;
    default: // avoid warning
      return 0.;
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

/*int ReactionFFABM_DISNC::Integrand_Cuhre(const int* ndim, const cubareal* x, const int *ncomp, cubareal* ff, void *userdata) {
  double f2(0), f2b(0), f2c(0), fl(0), flc(0), flb(0), f3(0), f3b(0), f3c(0);
  const integration_params& integrationParams = *(integration_params*)userdata;
  double xi = integrationParams.xi;
  double xip = xi + (1. - xi) * (x[0]);
  //printf("ReactionFFABM_DISNC::Integrand_Cuhre ndim = %d x[0] = %f ncomp = %d ff[0] = %f\n", *ndim, x[0], *ncomp, ff[0]);
  if(xip >= 1.) {
    ff[0] = 0.;
  }
  if (integrationParams.order == -1)
    sf_abkm_wrap_(xip, integrationParams.q2[integrationParams.i],
              f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b,
              integrationParams.ncflag, integrationParams.charge, integrationParams.polarity, *integrationParams._sin2thwPtr, integrationParams.cos2thw, *integrationParams._mzPtr);
  else
    sf_abkm_wrap_order_(xip, integrationParams.q2[integrationParams.i],
              f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b,
              integrationParams.ncflag, integrationParams.charge, integrationParams.polarity, *integrationParams._sin2thwPtr, integrationParams.cos2thw, *integrationParams._mzPtr, integrationParams.order);
  if (integrationParams.flag_calc_fl == 0) {
    if (integrationParams.flag_flavour == 1) {
      ff[0] = f2/xip/xip*(1.-xi);
      return 0;
    }
    else if (integrationParams.flag_flavour == 2) {
      ff[0] = f2c/xip/xip*(1.-xi);
      return 0;
    }
    else if (integrationParams.flag_flavour == 3) {
      ff[0] = f2b/xip/xip*(1.-xi);
      return 0;
    }
    else {
      ff[0] = 0.;
      return 0;
    }
  }
  else if (integrationParams.flag_calc_fl == 1) {
    if (integrationParams.flag_flavour == 1) {
      ff[0] = fl/xip/xip*(1.-xi);
      return 0;
    }
    else if (integrationParams.flag_flavour == 2) {
      ff[0] = flc/xip/xip*(1.-xi);
      return 0;
    }
    else if (integrationParams.flag_flavour == 3) {
      ff[0] = flb/xip/xip*(1.-xi);
      return 0;
    }
    else {
      ff[0] = 0.;
      return 0;
    }
  }
  else if (integrationParams.flag_calc_fl == 2) {
    if (integrationParams.flag_flavour == 1) {
      ff[0] = f3/xip/xip*(1.-xi);
      return 0;
    }
    else if (integrationParams.flag_flavour == 2) {
      ff[0] = f3c/xip/xip*(1.-xi);
      return 0;
    }
    else if (integrationParams.flag_flavour == 3) {
      ff[0] = f3b/xip/xip*(1.-xi);
      return 0;
    }
    else {
      ff[0] = 0.;
      return 0;
    }
  }
  throw 1;
};*/
