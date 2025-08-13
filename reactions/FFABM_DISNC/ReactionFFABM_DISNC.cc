
/*
   @file ReactionFFABM_DISNC.cc
   @date 2017-09-29
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-09-29
*/

#include "ReactionFFABM_DISNC.h"
#include "ABM.h"
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
double f2qcd_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2);
double flqcd_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2);
double f3qcd_(const int& nb, const int& nt, const int& ni, const double& xb, const double& q2);
double f2charm_ffn_(const double& xb, const double& q2, const int& nq);
double flcharm_ffn_(const double& xb, const double& q2, const int& nq);
// TODO to be removed
void sf_abkm_wrap_(const double &x, const double &q2,
  const double &f2abkm, const double &flabkm, const double &f3abkm,
  const double &f2cabkm, const double &flcabkm, const double &f3cabkm,
  const double &f2babkm, const double &flbabkm, const double &f3babkm,
  const int &ncflag, const double &charge, const double &polar,
  const double &sin2thw, const double &cos2thw, const double &MZ, const int &nt=1);
void abkm_set_input_(const int& kschemepdfin, const int& kordpdfin,
                     const double& rmass8in, const double& rmass10in, const int& msbarmin,
                     double& hqscale1in, const double& hqscale2in, const int& flord);
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
  _FLAG_FAST = td->hasParam("_FLAG_FAST") ? td->getParamI("_FLAG_FAST") : 0;

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

  printf("---------------------------------------------\n");
  printf("INFO from ReactionFFABM_DISNC:\n");
  printf("FF ABM running mass def? T(rue), (F)alse: %c\n", _msbarmin ? 'T' : 'F');
  printf("O(alpha_S) F_L - O(alpha_S) F2 = %d\n", _ordfl);
  printf("factorisation scale for heavy quarks  is set to sqrt(%f * Q^2 + %f * 4m_q^2\n", _hqscale1in, _hqscale2in);
  printf("---------------------------------------------\n");

  unsigned termID = td->id;
  _order[termID] = OrderMap(td->getParamS("Order")) - 1;
  _orderHQ[termID] = (td->hasParam("OrderHQ")) ? OrderMap(td->getParamS("OrderHQ")) - 1 : -1;
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

  abm::set_hq_masses(*_mcPtr, *_mbPtr);

  //if (_ht) {
  //  _ht->update();
  //}
  for (auto ht : _ht) {
    if (ht.second) {
      ht.second->update();
    }
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
    abm::pdffillgrid();

    // NC
    static constexpr int ncflag = 1;

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
        if (_FLAG_FAST) {
          if (GetDataFlav(dataSetID) == ReactionBaseDISNC::dataFlav::incl || GetDataFlav(dataSetID) == ReactionBaseDISNC::dataFlav::l) {
            //f2 = abm::calc_point_strfun(abm::SFproc::nc, abm::SFtype::f2, abm::SFflav::l, q2[i], x[i], -1, GetCharge(dataSetID), GetPolarisation(dataSetID), _sin2thetaW, *_mzPtr);
            f2 = calc_point_strfun(ReactionBaseDISNC::dataType::f2, ReactionBaseDISNC::dataFlav::l, q2[i], x[i], dataSetID, -1, GetCharge(dataSetID), _sin2thetaW, GetPolarisation(dataSetID), *_mzPtr);
            fl = calc_point_strfun(ReactionBaseDISNC::dataType::fl, ReactionBaseDISNC::dataFlav::l, q2[i], x[i], dataSetID, -1, GetCharge(dataSetID), _sin2thetaW, GetPolarisation(dataSetID), *_mzPtr);
            f3 = calc_point_strfun(ReactionBaseDISNC::dataType::f3, ReactionBaseDISNC::dataFlav::l, q2[i], x[i], dataSetID, -1, GetCharge(dataSetID), _sin2thetaW, GetPolarisation(dataSetID), *_mzPtr);
          }
          if (GetDataFlav(dataSetID) == ReactionBaseDISNC::dataFlav::incl || GetDataFlav(dataSetID) == ReactionBaseDISNC::dataFlav::c) {
            f2c = calc_point_strfun(ReactionBaseDISNC::dataType::f2, ReactionBaseDISNC::dataFlav::c, q2[i], x[i], dataSetID, -1, GetCharge(dataSetID), _sin2thetaW, GetPolarisation(dataSetID), *_mzPtr);
            flc = calc_point_strfun(ReactionBaseDISNC::dataType::fl, ReactionBaseDISNC::dataFlav::c, q2[i], x[i], dataSetID, -1, GetCharge(dataSetID), _sin2thetaW, GetPolarisation(dataSetID), *_mzPtr);
            f3c = calc_point_strfun(ReactionBaseDISNC::dataType::f3, ReactionBaseDISNC::dataFlav::c, q2[i], x[i], dataSetID, -1, GetCharge(dataSetID), _sin2thetaW, GetPolarisation(dataSetID), *_mzPtr);
          }
          if (GetDataFlav(dataSetID) == ReactionBaseDISNC::dataFlav::incl || GetDataFlav(dataSetID) == ReactionBaseDISNC::dataFlav::b) {
            f2b = calc_point_strfun(ReactionBaseDISNC::dataType::f2, ReactionBaseDISNC::dataFlav::b, q2[i], x[i], dataSetID, -1, GetCharge(dataSetID), _sin2thetaW, GetPolarisation(dataSetID), *_mzPtr);
            flb = calc_point_strfun(ReactionBaseDISNC::dataType::fl, ReactionBaseDISNC::dataFlav::b, q2[i], x[i], dataSetID, -1, GetCharge(dataSetID), _sin2thetaW, GetPolarisation(dataSetID), *_mzPtr);
            f3b = calc_point_strfun(ReactionBaseDISNC::dataType::f3, ReactionBaseDISNC::dataFlav::b, q2[i], x[i], dataSetID, -1, GetCharge(dataSetID), _sin2thetaW, GetPolarisation(dataSetID), *_mzPtr);
          }
        }
        else {
          abkm_set_input_(_kschemepdfin, _order[dataSetID], *_mcPtr, *_mbPtr, _msbarmin, _hqscale1in, _hqscale2in, _ordfl);
          sf_abkm_wrap_(x[i], q2[i], f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b, ncflag, charge, polarity, *_sin2thwPtr, cos2thw, *_mzPtr);
          //printf("%f %f %f\n", f2, f2c, f2b);
        }
                    
        double f3out_bar = 0.;
        if(_nuke[dataSetID] && _nuke[dataSetID]->need_f3bar()) {
          // need F3bar for nuclear corrections and antineutrino
          // we can calculate it now, because HT and TMC (calculated later) do not apply to F3
          int charge_bar = -1 * charge;
          double f2_bar(0), f2b_bar(0), f2c_bar(0), fl_bar(0), flc_bar(0), flb_bar(0), f3_bar(0), f3c_bar(0), f3b_bar(0);
          if (_FLAG_FAST) {
            if (GetDataFlav(dataSetID) == ReactionBaseDISNC::dataFlav::incl || GetDataFlav(dataSetID) == ReactionBaseDISNC::dataFlav::l) {
              f3_bar = calc_point_strfun(ReactionBaseDISNC::dataType::f3, ReactionBaseDISNC::dataFlav::l, q2[i], x[i], dataSetID, -1, charge_bar, _sin2thetaW, GetPolarisation(dataSetID), *_mzPtr);
            }
            if (GetDataFlav(dataSetID) == ReactionBaseDISNC::dataFlav::incl || GetDataFlav(dataSetID) == ReactionBaseDISNC::dataFlav::c) {
              f3c_bar = calc_point_strfun(ReactionBaseDISNC::dataType::f3, ReactionBaseDISNC::dataFlav::c, q2[i], x[i], dataSetID, -1, charge_bar, _sin2thetaW, GetPolarisation(dataSetID), *_mzPtr);
            }
            if (GetDataFlav(dataSetID) == ReactionBaseDISNC::dataFlav::incl || GetDataFlav(dataSetID) == ReactionBaseDISNC::dataFlav::b) {
              f3b_bar = calc_point_strfun(ReactionBaseDISNC::dataType::f3, ReactionBaseDISNC::dataFlav::b, q2[i], x[i], dataSetID, -1, charge_bar, _sin2thetaW, GetPolarisation(dataSetID), *_mzPtr);
            }
          }
          else {
            sf_abkm_wrap_(x[i], q2[i], f2_bar, fl_bar, f3_bar, f2c_bar, flc_bar, f3c_bar, f2b_bar, flb_bar, f3b_bar, ncflag, charge_bar, polarity, *_sin2thwPtr, cos2thw, *_mzPtr);
          }
          //double f3_bar(0), f3c_bar(0), f3b_bar(0);
          f3out_bar = x[i] * combine_flavours(GetDataFlav(dataSetID), f3_bar, f3c_bar, f3b_bar);
        }
        if (_tmc[dataSetID]) {
          const bool flag_fl = true;
          const bool flag_f3 = false;
          //const bool flag_f3 = true; [not implemented]
          _tmc[dataSetID]->apply(f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b, flag_fl, flag_f3, q2[i], x[i], ncflag, charge, polarity, cos2thw, *_mzPtr, this);
        }
        //if(_flag_ht[dataSetID]) {
        //  _ht->apply(q2[i], x[i], f2, fl);
        //}      
        if(_ht[dataSetID]) {
          _ht[dataSetID]->apply(q2[i], x[i], f2, fl);
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

double ReactionFFABM_DISNC::calc_point_strfun(const ReactionBaseDISNC::dataType ftype, const ReactionBaseDISNC::dataFlav flav, const double q2, const double x, const int dataSetID, const int order, const int charge, const double sin2thw, const double polar, const double mz) {
  int orderALL = (order >= 0) ? order : _order[dataSetID];
  int orderHQ = (order >= 0) ? order : _orderHQ[dataSetID];
  int orderFL = (order >= 0) ? order : _ordfl;
  abm::set_scheme_and_order(_kschemepdfin, orderALL, _msbarmin, orderFL, orderHQ);
  static constexpr int nt = 1; // proton
  switch (flav) {
    case ReactionBaseDISNC::dataFlav::l: {
      static constexpr double eleAxial = -0.5;
      const double eleVec = -0.5 + 2 * sin2thw;
      const double facgz = - eleVec - charge * polar * eleAxial;
      const double faczz = eleVec * eleVec + eleAxial * eleAxial + 2 * charge * polar * eleAxial * eleVec;
      const double facgzf3 = -1 * eleAxial - charge * polar * eleVec;
      const double faczzf3 = 2 * eleAxial * eleVec + charge * polar * (eleVec * eleVec + eleAxial * eleAxial);
      const double PZ = 1. / (4 * sin2thw * (1 - sin2thw) * (1 + mz * mz / q2));
      switch (ftype) {
        case ReactionBaseDISNC::dataType::f2:
          return f2qcd_(3, nt, 22, x, q2) + facgz * PZ * f2qcd_(3, nt, 25, x, q2) + faczz * PZ * PZ * f2qcd_(3, nt, 23, x, q2);
        case ReactionBaseDISNC::dataType::fl:
          return flqcd_(3, nt, 22, x, q2) + facgz * PZ * flqcd_(3, nt, 25, x, q2) + faczz * PZ * PZ * flqcd_(3, nt, 23, x, q2);
        case ReactionBaseDISNC::dataType::f3:
          return -1 * charge * (facgzf3 * PZ * f3qcd_(3, nt, 25, x, q2) + faczzf3 * PZ * PZ * f3qcd_(3, nt, 23, x, q2));
      }
    }
    case ReactionBaseDISNC::dataFlav::c:
      switch (ftype) {
        case ReactionBaseDISNC::dataType::f2:
          return f2charm_ffn_(x, q2, 8);
        case ReactionBaseDISNC::dataType::fl:
          return flcharm_ffn_(x, q2, 8);
        case ReactionBaseDISNC::dataType::f3:
          return 0.;
      }
    case ReactionBaseDISNC::dataFlav::b:
      switch (ftype) {
        case ReactionBaseDISNC::dataType::f2:
          return f2charm_ffn_(x, q2, 10);
        case ReactionBaseDISNC::dataType::fl:
          return flcharm_ffn_(x, q2, 10);
        case ReactionBaseDISNC::dataType::f3:
          return 0.;
      }
  }
  hf_errlog(28022501, "F: Unsupported structure function type or flavour");
  return 0; // avoid warning
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
    default:
      hf_errlog(28022501, "F: Unsupported flavour");
      return 0.; // avoid warning
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
