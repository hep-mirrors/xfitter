
#include "ReactionBaseFFABM.h"
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
#include <iterator>
#include <algorithm>


// Initialize at the start of the computation
//void ReactionBaseFFABM::atStart()
//{
//  printf("ReactionBaseFFABM::atStart()\n");
  // do not call parent atStart(): it initialises QCDNUM
  // Super::atStart();
//  ReactionTheory::atStart();
//}

void ReactionBaseFFABM::initTerm(TermData *td)
{
  printf("ReactionBaseFFABM::initTerm()\n");
  Super::initTerm(td);
  unsigned termID = td->id;

  // scales mu^2 = scalea1 * Q^2 + scaleb1 * 4*m_h^2 (default scalea1 = scaleb1 = 1.0)
  const double hqscale1in = *td->getParamD("scalea1");
  const double hqscale2in = *td->getParamD("scaleb1");
  abm::set_hq_scales(hqscale1in, hqscale2in);

  // pole or MCbar running mass treatment (default pole)
  const int msbarmin = td->getParamI("runm");
  // O(alpha_S) F_L = O(alpha_S) F_2 + ordfl (default ordfl = 1)
  const int ordfl = td->getParamI("ordfl");
  const int order = OrderMap(td->getParamS("Order")) - 1;
  const auto& str_orderhq = td->getParamS("OrderHQ");
  const int orderHQ = str_orderhq == "" ? -1 : OrderMap(str_orderhq) - 1;

  // control x range (certain PDF sets have limited x_min, x_max)
  abm::set_xbmin(*td->getParamD("xbmin"));
  abm::set_xbmax(*td->getParamD("xbmax"));
  _q2mincomp = *td->getParamD("q2mincomp");
  abm::initgridconst();

  // heavy quark masses
  _mcPtr = td->getParamD("mch");
  _mbPtr = td->getParamD("mbt");

  printf("---------------------------------------------\n");
  printf("INFO from %s: ", getReactionName().c_str());
  printf("running mass = %c", msbarmin ? 'T' : 'F');
  printf(", order = %d", order);
  printf(", order HQ = %d", orderHQ);
  printf(", O(alpha_S) F_L - O(alpha_S) F2 = %d", ordfl);
  printf(", factorisation scale for heavy quarks = sqrt(%f * Q^2 + %f * 4m_q^2\n", hqscale1in, hqscale2in);
  printf("---------------------------------------------\n");

  auto nBins = td->getNbins();
  //if(_integrated.find(termID) != _integrated.end())
  //  nBins = _integrated[termID]->getBinValuesQ2()->size();
  _f2abm[termID].resize(nBins);
  _flabm[termID].resize(nBins);
  _f3abm[termID].resize(nBins);

  const double* mzPtr = td->getParamD("Mz");
  const double* sin2thwPtr = td->getParamD("sin2thW");

  const auto& q2 = *GetBinValues(td, "Q2");
  const auto& x = *GetBinValues(td, "x");
  // groups for parallel
  const auto& parallel_task_distribution = td->getParamS("parallel_task_distribution");
  _task_distr = ForkPool::TaskDistribution::chunky;
  if(parallel_task_distribution == "chunky") {
    _task_distr = ForkPool::TaskDistribution::chunky;
  }
  else if(parallel_task_distribution == "cyclic") {
    _task_distr = ForkPool::TaskDistribution::cyclic;
  }
  else {
    auto msg = "F: unknown task distribution " + parallel_task_distribution;
    hf_errlog(2026061601,msg);
  }
  const auto& group = td->getParamS("group_parallel");
  printf("group = %s\n", group.c_str());
  size_t offset = 0;
  size_t nBinsActive = std::count_if(std::begin(q2), std::end(q2), [this](double q2) { return q2 >= this->_q2mincomp; });
  auto it = _grouped_data_points.find(group);
  if(it == _grouped_data_points.end()) {
    printf("inserting group = %s with nBinsActive = %lu\n", group.c_str(), nBinsActive);
    it = _grouped_data_points.insert(std::make_pair(group, std::vector<DataPoint>(nBinsActive))).first;
  }
  else {
    size_t size0 = it->second.size();
    printf("found group = %s size0,nBinsActive = %lu,%lu\n", group.c_str(), size0, nBinsActive);
    it->second.resize(size0 + nBinsActive);
    offset = size0;
  }
  auto& vec2 = it->second;
  int iActive = 0;
  for(int i = 0; i < nBins; i++) {
    if(q2[i] <= _q2mincomp) {
      _f2abm[i] = 0.0;
      _flabm[i] = 0.0;
      _f3abm[i] = 0.0;
      continue;
    }
    //vec2[offset + i] = vec[i];
    auto& point = vec2[offset + iActive];
    point.td = td;
    point.datasetID = td->id;
    point.i = i;
    point.proc_NCCC = GetProc();
    point.flav = dataFlav(GetDataFlav(td->id));
    point.ord = order;
    point.ordHQ = orderHQ;
    point.ordFL = order;
    point.msbarmin = msbarmin;
    point.charge = GetCharge(td->id);
    point.polar = GetPolarisation(td->id);
    point.sin2thetaWPtr = sin2thwPtr;
    point.mz = mzPtr;
    point.q2 = q2[i];
    point.x = x[i];
    point.ht = getHT(td->id);
    point.tmc = getTMC(td->id);
    point.nuke = getNUKE(td->id);
    point.f2 = 0.0;
    point.fl = 0.0;
    point.f3 = 0.0;
    iActive++;
  }
}

//
void ReactionBaseFFABM::atIteration()
{
  printf("ReactionBaseFFABM::atIteration()\n");
  Super::atIteration();
  abm::set_hq_masses(*_mcPtr, *_mbPtr);
  _grouped_data_points.begin()->second[0].td->actualizeWrappers();
  abm::pdffillgrid();

  // parallel computaion for groups of data points
  printf("atIteration ReactionTheory::_ncpu = %d _flagComputeAtIteration = %d\n", ReactionTheory::_ncpu, _flagComputeAtIteration);
  if(ReactionTheory::_ncpu == 1) {
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
      bool need_pdffillgrid = vec[0].td->actualizeWrappers();
      printf("checking pdffillgrid() group = %s, datasetID = %d\n", it.first.c_str(), vec[0].datasetID);
      if(need_pdffillgrid) {
        printf("extra pdffillgrid() group = %s, datasetID = %d\n", it.first.c_str(), vec[0].datasetID);
        abm::pdffillgrid();
      }
      ForkPool pool(ReactionTheory::_ncpu, _task_distr);
      int np = vec.size();
      ForkPool::SharedMemory shm(sizeof(double) * 3 * np + sizeof(int) * np * 2);
      double* f2 = shm.data<double>();
      double* fl = f2 + np;
      double* f3 = fl + np;
      int* datasetID = (int*)(f3 + np);
      int* i_orig = datasetID + np;
      pool.parallel_for(vec, [&vec, f2, fl, f3, datasetID, i_orig](size_t i) {
        auto& point = vec[i];
        point.calc();
        f2[i] = point.f2;
        fl[i] = point.fl;
        f3[i] = point.f3;
        datasetID[i] = point.datasetID;
        i_orig[i] = point.i;
      });
      for(size_t i = 0; i < np; i++) {
        _f2abm[datasetID[i]][i_orig[i]] = f2[i];
        _flabm[datasetID[i]][i_orig[i]] = fl[i];
        _f3abm[datasetID[i]][i_orig[i]] = f3[i];
      }
    }
  }
}

// Calculates one data point at (Q2,x) and returns values f2, fl, f3
void ReactionBaseFFABM::DataPoint::calc()
{
  auto combine_flavours = [](const dataFlav flav, const double f, const double fc, const double fb) {
    switch (flav)
    {
      case dataFlav::incl:
        return f + fc + fb;
      case dataFlav::c:
        return fc;
      case dataFlav::b:
        return fb;
      default:
        hf_errlog(28022501, "F: Unsupported flavour");
        return 0.; // avoid warning
    }
  };

  bool need_pdffillgrid = td->actualizeWrappers();
  if(need_pdffillgrid) {
    printf("extra pdffillgrid() by datasetID = %d point = %d (proc_NCCC = %d)\n", datasetID, i, (int)proc_NCCC);
    //fflush(stdout);
    abm::pdffillgrid();
  }
  //printf("q2,x = %f,%f\n", q2, x);fflush(stdout);
  f2 = fl = f3 = 0.;
  double f2l(0), f2b(0), f2c(0), fll(0), flc(0), flb(0), f3l(0), f3b(0), f3c(0);
  if (flav == dataFlav::incl || flav == dataFlav::l) {
    f2l = abm::calc_point_strfun(proc_NCCC, abm::SFtype::f2, abm::SFflav::l, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, *mz);
    fll = abm::calc_point_strfun(proc_NCCC, abm::SFtype::fl, abm::SFflav::l, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, *mz);
    f3l = abm::calc_point_strfun(proc_NCCC, abm::SFtype::f3, abm::SFflav::l, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, *mz);
  }
  if (flav == dataFlav::incl || flav == dataFlav::c) {
    f2c = abm::calc_point_strfun(proc_NCCC, abm::SFtype::f2, abm::SFflav::c, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, *mz);
    flc = abm::calc_point_strfun(proc_NCCC, abm::SFtype::fl, abm::SFflav::c, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, *mz);
    f3c = abm::calc_point_strfun(proc_NCCC, abm::SFtype::f3, abm::SFflav::c, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, *mz);
  }
  if (flav == dataFlav::incl || flav == dataFlav::b) {
    f2b = abm::calc_point_strfun(proc_NCCC, abm::SFtype::f2, abm::SFflav::b, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, *mz);
    flb = abm::calc_point_strfun(proc_NCCC, abm::SFtype::fl, abm::SFflav::b, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, *mz);
    f3b = abm::calc_point_strfun(proc_NCCC, abm::SFtype::f3, abm::SFflav::b, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, *mz);
  }
                
  double f3lout_bar = 0.;
  if(nuke && nuke->need_f3bar()) {
    // need F3bar for nuclear corrections and antineutrino
    // we can calculate it now, because HT and TMC (calculated later) do not apply to F3
    int charge_bar = -1 * charge;
    double f2l_bar(0), f2b_bar(0), f2c_bar(0), fll_bar(0), flc_bar(0), flb_bar(0), f3l_bar(0), f3c_bar(0), f3b_bar(0);
    if (flav == dataFlav::incl || flav == dataFlav::l) {
      f3l_bar = abm::calc_point_strfun(proc_NCCC, abm::SFtype::f3, abm::SFflav::l, q2, x, -1, ord, ordHQ, ordHQ, msbarmin, charge_bar, *sin2thetaWPtr, polar, *mz);
    }
    if (flav == dataFlav::incl || flav == dataFlav::c) {
      f3c_bar = abm::calc_point_strfun(proc_NCCC, abm::SFtype::f3, abm::SFflav::c, q2, x, -1, ord, ordHQ, ordHQ, msbarmin, charge_bar, *sin2thetaWPtr, polar, *mz);
    }
    if (flav == dataFlav::incl || flav == dataFlav::b) {
      f3b_bar = abm::calc_point_strfun(proc_NCCC, abm::SFtype::f3, abm::SFflav::b, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge_bar, *sin2thetaWPtr, polar, *mz);
    }
    f3lout_bar = x * combine_flavours(flav, f3l_bar, f3c_bar, f3b_bar);
  }
  if (tmc) {
    static constexpr abm::SFproc ncflag = abm::SFproc::nc;
    tmc->apply(f2l, fll, f3l, f2c, flc, f3c, f2b, flb, f3b, q2, x, ncflag, ord, ordHQ, ordFL, msbarmin, charge, polar, 1.-*sin2thetaWPtr, *mz);
  }
  if(ht) {
    // HT is applied only to F2 and FL light flavour part
    ht->apply(q2, x, f2l, fll, f3l);
  }
  f2 = combine_flavours(flav, f2l, f2c, f2b);
  fl = combine_flavours(flav, fll, flc, flb);
  f3 = x * combine_flavours(flav, f3l, f3c, f3b);
  // apply nuclear corrections to the sum of light+c+b because corrections for charm and non-charm (kint=4,5) are not implemented
  if (nuke) {
    nuke->apply(q2, x, f2, fl, f3);
  }
}
