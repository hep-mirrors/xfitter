
/*
   @file ReactionFFABM_DISNC.cc
   @date 2017-09-29
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-09-29
*/

#include "ReactionFFABM_DISNC.h"
#include "DIS_HT.h"
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

    // Calculate i-th data point, return values f2, fl, f3
    auto calc_point = [&](int i, double& f2out, double& flout, double& f3out) {
      f2out = 0.;
      flout = 0.;
      f3out = 0.;
      if (q2[i] > 1.0)
      {
        //printf("q2,x = %f,%f\n", q2[i], x[i]);
        sf_abkm_wrap_(x[i], q2[i],
                      f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b,
                      ncflag, charge, polarity, *_sin2thwPtr, cos2thw, *_mzPtr);
        if(_flag_tmc[dataSetID]) {
          // TODO ???
          if ((_tmc_xmin[i] == 0. || _tmc_xmin[i] < x[i]) && (_tmc_logxlogq2min[i] == 0. || _tmc_logxlogq2min[i] < log(x[i])*log(q2[i])))
            apply_tmc(_tmc_integration_method[dataSetID], f2, fl, f3, 1, q2, x, ncflag, charge, polarity, cos2thw, i);
        }
        if(_flag_ht[dataSetID]) {
          _ht->apply(q2[i], x[i], f2, fl);
        }      
      }
      switch (GetDataFlav(dataSetID))
      {
        case dataFlav::incl:
          f2out = f2 + f2c + f2b;
          flout = fl + flc + flb;
          f3out = x[i] * (f3 + f3c + f3b);
          break;
        case dataFlav::c:
          f2out = f2c;
          flout = flc;
          f3out = x[i] * f3c;
          break;
        case dataFlav::b:
          f2out = f2b;
          flout = flb;
          f3out = x[i] * f3b;
          break;
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

double ReactionFFABM_DISNC::apply_tmc(const int method, double& f2, double& fl, double& f3, const int flag_flavour, const std::valarray<double>& q2, const std::valarray<double>& x,
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
