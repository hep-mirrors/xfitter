
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
#include "ForkPool.h"
#include <numeric>


// the class factories
extern "C" ReactionFFABM_DISNC *create()
{
  return new ReactionFFABM_DISNC();
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

  printf("---------------------------------------------\n");
  printf("INFO from ReactionFFABM_DISNC:\n");
  printf("FF ABM running mass def? T(rue), (F)alse: %c\n", _msbarmin[termID] ? 'T' : 'F');
  printf("Order = %d\n", _order[termID]);
  printf("Order HQ = %d\n", _orderHQ[termID]);
  printf("O(alpha_S) F_L - O(alpha_S) F2 = %d\n", _ordfl[termID]);
  printf("factorisation scale for heavy quarks  is set to sqrt(%f * Q^2 + %f * 4m_q^2\n", _hqscale1in, _hqscale2in);
  printf("---------------------------------------------\n");

  auto nBins = td->getNbins();
  if(_integrated.find(termID) != _integrated.end())
    nBins = _integrated[termID]->getBinValuesQ2()->size();
  _f2abm[termID].resize(nBins);
  _flabm[termID].resize(nBins);
  _f3abm[termID].resize(nBins);

  _mzPtr = td->getParamD("Mz");
  _sin2thwPtr = td->getParamD("sin2thW");

  // parallel
  _ncpu[td->id] = getNCPU(td);

  _data_points[td->id] = std::vector<DataPoint>(nBins);
  auto& vec = _data_points[td->id];
  const auto& q2 = *GetBinValues(td, "Q2");
  const auto& x = *GetBinValues(td, "x");
  for(int i = 0; i < nBins; i++) {
    auto& point = vec[i];
    point.datasetID = td->id;
    point.i = i;
    point.flav = GetDataFlav(td->id);
    point.ord = _order[td->id];
    point.ordHQ = _orderHQ[td->id];
    point.ordFL = _order[td->id];
    point.msbarmin = _msbarmin[td->id];
    point.charge = GetCharge(td->id);
    point.charge = GetPolarisation(td->id);
    point.sin2thetaW = _sin2thetaW;
    point.mz = *_mzPtr;
    point.q2 = q2[i];
    point.x = x[i];
    point.f2 = 0.0;
    point.fl = 0.0;
    point.f3 = 0.0;
  }
  // groups for parallel
  std::string group = "default";
  if(td->hasParam("group_parallel")) {
    group = td->getParamS("group_parallel");
  }
  auto it = _grouped_data_points.find(group);
  size_t offset = 0;
  if(it == _grouped_data_points.end()) {
    printf("nBins = %lu\n", nBins);fflush(stdout);
    it = _grouped_data_points.insert(std::make_pair(group, std::vector<DataPoint>(nBins))).first;
  }
  else {
    size_t size0 = it->second.size();
    printf("size0,nBins = %lu,%lu\n", size0, nBins);fflush(stdout);
    it->second.resize(size0 + nBins);
    offset = size0;
  }
  auto& vec2 = it->second;
  for(int i = 0; i < nBins; i++) {
    vec2[offset + i] = vec[i];
  }
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

  auto td = _tdDS.begin()->second;
  //td->actualizeWrappers();
  //abm::pdffillgrid();
  _need_pdffillgrid = true;

  // parallel computaion for groups of data points
  printf("atIteration\n");fflush(stdout);
  if (ReactionTheory::_ncpu > 1 || _flagComputeAtIteration) {
    if(ReactionTheory::_ncpu) {
      for(auto& it : _grouped_data_points) {
        for(auto& point : it.second) {
          point.calc();
          _f2abm[point.datasetID][point.i] = point.f2;
          _flabm[point.datasetID][point.i] = point.fl;
          _f3abm[point.datasetID][point.i] = point.f3;
        }
      }
    }
    else {
      for(auto& it : _grouped_data_points) {
        auto& vec = it.second;
        ForkPool pool(ReactionTheory::_ncpu);
        int np = vec.size();
        ForkPool::SharedMemory shm(sizeof(double) * 3 * np);
        double* f2 = shm.data<double>();
        double* fl = f2 + np;
        double* f3 = fl + np;
        pool.parallel_for(vec, [&vec](size_t i) {vec[i].calc();});
        for (auto& point : vec) {
          _f2abm[point.datasetID][point.i] = point.f2;
          _flabm[point.datasetID][point.i] = point.fl;
          _f3abm[point.datasetID][point.i] = point.f3;
        }
      }
    }
  }
}

// Place calculations in one function, to optimize calls.
void ReactionFFABM_DISNC::calcF2FL(unsigned dataSetID)
{
  printf("calcF2FL\n");fflush(stdout);
  if ((_f2abm[dataSetID][0] < -99.))
  { // compute
    auto td = _tdDS[dataSetID];
    bool need_update = td->actualizeWrappers();
    if (_need_pdffillgrid || need_update) {
      printf("ReactionFFABM_DISNC pdffillgrid()\n");
      abm::pdffillgrid();
      _need_pdffillgrid = false;
    }

    double charge = GetCharge(dataSetID);
    double polarity = GetPolarisation(dataSetID);

    // Get x,Q2 arrays:
    auto *q2p = GetBinValues(td, "Q2"), *xp = GetBinValues(td, "x");
    auto q2 = *q2p, x = *xp;

    // Number of data points
    // SZ perhaps GetNpoint does not work for integrated cross sections (not tested)
    //const size_t Np = GetNpoint(dataSetID);
    const size_t np = xp->size();

    double f2(0), f2b(0), f2c(0), fl(0), flc(0), flb(0), f3(0), f3b(0), f3c(0);
    //double cos2thw = 1.0 - *_sin2thwPtr;

    //printf("FFABM _ncpu = %d\n", _ncpu[dataSetID]);
    int ncpu =  xfitter::xf_ncpu(_ncpu[dataSetID]);
    //printf("FFABM ncpu = %d\n", ncpu);

    //auto calc_point = [&](int i, double& f2out, double& flout, double& f3out) {
    //  _data_points[dataSetID][i].calc();
    //  f2out = f2;
    //  flout = fl;
    //  f3out = f3;
    //};

    auto& vec = _data_points[dataSetID];
    if (ncpu == 1) {
      for (size_t i = 0; i < np; i++) {
        auto& point = vec[i];
        point.calc();
        _f2abm[dataSetID][i] = point.f2;
        _flabm[dataSetID][i] = point.fl;
        _f3abm[dataSetID][i] = point.f3;
      }
    }
    else {
      if(1) {ForkPool pool(ncpu);
      ForkPool::SharedMemory shm(sizeof(double) * 3 * np);
      double* f2 = shm.data<double>();
      double* fl = f2 + np;
      double* f3 = fl + np;
      pool.parallel_for(vec, [&vec](size_t i) {vec[i].calc();});
      for (auto& v : vec) {
        _f2abm[dataSetID][v.i] = v.f2;
        _flabm[dataSetID][v.i] = v.fl;
        _f3abm[dataSetID][v.i] = v.f3;
      }}
      /*if(0) {// Shared memory for predictions
      int shmid;
      double* sharedArray;
      shmid = shmget(IPC_PRIVATE, sizeof(double) * np * 3, IPC_CREAT | 0666);
      if (shmid < 0) {
        hf_errlog(2023060200,"F: Failed to create shared memory segment");
      }
      sharedArray = static_cast<double*>(shmat(shmid, nullptr, 0));
      if (sharedArray == reinterpret_cast<double*>(-1)) {
        hf_errlog(2023060201,"F: Failed to attach shared memory segment");
      }
      // define Chunks
      int chunkSize = np / ncpu;
      int reminder  = np % ncpu; 
      int first = 0;
      int startIndex = 0;
      int endIndex = 0;
      // loop over all
      for (int icpu = 0; icpu < std::min(ncpu, int(np)); icpu++) {
        startIndex = endIndex;
        endIndex   = startIndex + chunkSize;
        if (icpu < reminder) {
          endIndex += 1;
        }
        printf("icpu = %d startIndex = %d endIndex = %d\n", icpu, startIndex, endIndex);
        pid_t pid = xfitter::xf_fork( std::min(ncpu, int(np))  );
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
            sharedArray[i+np] = fl;
            sharedArray[i+2*np] = f3;
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
      for (size_t i = 0; i<np; i++) {
        _f2abm[dataSetID][i] = sharedArray[i];
        _flabm[dataSetID][i] = sharedArray[i+np];
        _f3abm[dataSetID][i] = sharedArray[i+2*np];
      }    
      // Detach and remove shared memory segments
      shmdt(sharedArray);
      shmctl(shmid, IPC_RMID, NULL);}*/
    }
  }
}

// Calculates one data point at (Q2,x) and returns values f2, fl, f3
void ReactionFFABM_DISNC::DataPoint::calc()
{
  f2 = fl = f3 = 0.;
  double f2l(0), f2b(0), f2c(0), fll(0), flc(0), flb(0), f3l(0), f3b(0), f3c(0);
  if (flav == ReactionBaseDISNC::dataFlav::incl || flav == ReactionBaseDISNC::dataFlav::l) {
    f2 = abm::calc_point_strfun_NC(abm::SFtype::f2, abm::SFflav::l, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, sin2thetaW, polar, mz);
    fl = abm::calc_point_strfun_NC(abm::SFtype::fl, abm::SFflav::l, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, sin2thetaW, polar, mz);
    f3 = abm::calc_point_strfun_NC(abm::SFtype::f3, abm::SFflav::l, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, sin2thetaW, polar, mz);
  }
  if (flav == ReactionBaseDISNC::dataFlav::incl || flav == ReactionBaseDISNC::dataFlav::c) {
    f2c = abm::calc_point_strfun_NC(abm::SFtype::f2, abm::SFflav::c, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, sin2thetaW, polar, mz);
    flc = abm::calc_point_strfun_NC(abm::SFtype::fl, abm::SFflav::c, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, sin2thetaW, polar, mz);
    f3c = abm::calc_point_strfun_NC(abm::SFtype::f3, abm::SFflav::c, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, sin2thetaW, polar, mz);
  }
  if (flav == ReactionBaseDISNC::dataFlav::incl || flav == ReactionBaseDISNC::dataFlav::b) {
    f2b = abm::calc_point_strfun_NC(abm::SFtype::f2, abm::SFflav::b, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, sin2thetaW, polar, mz);
    flb = abm::calc_point_strfun_NC(abm::SFtype::fl, abm::SFflav::b, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, sin2thetaW, polar, mz);
    f3b = abm::calc_point_strfun_NC(abm::SFtype::f3, abm::SFflav::b, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, sin2thetaW, polar, mz);
  }
                
  double f3lout_bar = 0.;
  if(nuke && nuke->need_f3bar()) {
    // need F3bar for nuclear corrections and antineutrino
    // we can calculate it now, because HT and TMC (calculated later) do not apply to F3
    int charge_bar = -1 * charge;
    double f2l_bar(0), f2b_bar(0), f2c_bar(0), fll_bar(0), flc_bar(0), flb_bar(0), f3l_bar(0), f3c_bar(0), f3b_bar(0);
    if (flav == ReactionBaseDISNC::dataFlav::incl || flav == ReactionBaseDISNC::dataFlav::l) {
      f3l_bar = abm::calc_point_strfun_NC(abm::SFtype::f3, abm::SFflav::l, q2, x, -1, ord, ordHQ, ordHQ, msbarmin, charge_bar, sin2thetaW, polar, mz);
    }
    if (flav == ReactionBaseDISNC::dataFlav::incl || flav == ReactionBaseDISNC::dataFlav::c) {
      f3c_bar = abm::calc_point_strfun_NC(abm::SFtype::f3, abm::SFflav::c, q2, x, -1, ord, ordHQ, ordHQ, msbarmin, charge_bar, sin2thetaW, polar, mz);
    }
    if (flav == ReactionBaseDISNC::dataFlav::incl || flav == ReactionBaseDISNC::dataFlav::b) {
      f3b_bar = abm::calc_point_strfun_NC(abm::SFtype::f3, abm::SFflav::b, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge_bar, sin2thetaW, polar, mz);
    }
    f3lout_bar = x * ReactionFFABM_DISNC::combine_flavours(flav, f3l_bar, f3c_bar, f3b_bar);
  }
  if (tmc) {
    static constexpr abm::SFproc ncflag = abm::SFproc::nc;
    tmc->apply(f2, fl, f3, f2c, flc, f3c, f2b, flb, f3b, q2, x, ncflag, ord, ordHQ, ordFL, msbarmin, charge, polar, 1.-sin2thetaW, mz);
  }
  if(ht) {
    // HT is applied only to F2 and FL light flavour part
    ht->apply(q2, x, f2, fl, f3);
  }      
  f2 = ReactionFFABM_DISNC::combine_flavours(flav, f2, f2c, f2b);
  fl = ReactionFFABM_DISNC::combine_flavours(flav, fl, flc, flb);
  f3 = x * combine_flavours(flav, f3, f3c, f3b);
  // apply nuclear corrections to the sum of light+c+b because corrections for charm and non-charm (kint=4,5) are not implemented
  if (nuke) {
    nuke->apply(q2, x, f2, fl, f3);
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
