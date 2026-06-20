
#include "ReactionFFABM_DISNC_CC.h"
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
#include <iterator>
#include <algorithm>


// the class factories
extern "C" ReactionFFABM_DISNC_CC *create()
{
  return new ReactionFFABM_DISNC_CC();
}

// Initialize at the start of the computation
void ReactionFFABM_DISNC_CC::atStart()
{
  // do not call parent atStart(): it initialises QCDNUM
  // Super::atStart();
  ReactionTheory::atStart();
}

void ReactionFFABM_DISNC_CC::initTerm(TermData *td)
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
  if(td->hasParam("q2mincomp"))
    _q2mincomp = *td->getParamD("q2mincomp");

  abm::initgridconst();

  // heavy quark masses
  _mcPtr = td->getParamD("mch");
  _mbPtr = td->getParamD("mbt");

  printf("---------------------------------------------\n");
  printf("INFO from ReactionFFABM_DISNC_CC:\n");
  printf("FF ABM running mass def? T(rue), (F)alse: %c\n", _msbarmin[termID] ? 'T' : 'F');
  printf("Order = %d\n", _order[termID]);
  printf("Order HQ = %d\n", _orderHQ[termID]);
  printf("O(alpha_S) F_L - O(alpha_S) F2 = %d\n", _ordfl[termID]);
  printf("factorisation scale for heavy quarks  is set to sqrt(%f * Q^2 + %f * 4m_q^2\n", _hqscale1in, _hqscale2in);
  printf("---------------------------------------------\n");

  auto nBins = td->getNbins();
  //if(_integrated.find(termID) != _integrated.end())
  //  nBins = _integrated[termID]->getBinValuesQ2()->size();
  _f2abm[termID].resize(nBins);
  _flabm[termID].resize(nBins);
  _f3abm[termID].resize(nBins);

  _mzPtr = td->getParamD("Mz");
  _sin2thwPtr = td->getParamD("sin2thW");

  // parallel
  _ncpu[td->id] = getNCPU(td);

  const auto& q2 = *GetBinValues(td, "Q2");
  const auto& x = *GetBinValues(td, "x");
  size_t nBinsActive = std::count_if(std::begin(q2), std::end(q2), [this](double q2) { return q2 >= this->_q2mincomp; });
  _data_points[td->id] = std::vector<DataPoint>(nBinsActive);
  auto& vec = _data_points[td->id];
  size_t iActive = 0;
  for(int i = 0; i < nBins; i++) {
    if(q2[i] <= _q2mincomp) {
      _f2abm[i] = 0.0;
      _flabm[i] = 0.0;
      _f3abm[i] = 0.0;
      continue;
    }
    auto& point = vec[iActive];
    point.datasetID = td->id;
    point.i = i;
    point.flav = dataFlav(GetDataFlav(td->id));
    point.ord = _order[td->id];
    point.ordHQ = _orderHQ[td->id];
    point.ordFL = _order[td->id];
    point.msbarmin = _msbarmin[td->id];
    point.charge = GetCharge(td->id);
    point.polar = GetPolarisation(td->id);
    point.sin2thetaWPtr = _sin2thwPtr;
    point.mz = *_mzPtr;
    point.q2 = q2[i];
    point.x = x[i];
    point.ht = _ht[td->id];
    point.tmc = _tmc[td->id];
    point.nuke = _nuke[td->id];
    point.f2 = 0.0;
    point.fl = 0.0;
    point.f3 = 0.0;
    iActive++;
  }
  // groups for parallel
  _task_distr = ForkPool::TaskDistribution::chunky;
  if(td->hasParam("parallel_task_distribution")) {
    auto str = td->getParamS("parallel_task_distribution");
    if(str == "chunky") {
      _task_distr = ForkPool::TaskDistribution::chunky;
    }
    else if(str == "cyclic") {
      _task_distr = ForkPool::TaskDistribution::cyclic;
    }
    else {
      auto msg = "F: unknown task distribution " + str;
      hf_errlog(2026061601,msg);
    }
  }
  std::string group = "default";
  if(td->hasParam("group_parallel")) {
    group = td->getParamS("group_parallel");
  }
  size_t offset = 0;
  auto it = _grouped_data_points.find(group);
  if(it == _grouped_data_points.end()) {
    printf("inserting group = %s with nBinsActive = %lu\n", group.c_str(), nBinsActive);fflush(stdout);
    it = _grouped_data_points.insert(std::make_pair(group, std::vector<DataPoint>(nBinsActive))).first;
  }
  else {
    size_t size0 = it->second.size();
    printf("found group = %s size0,nBinsActive = %lu,%lu\n", group.c_str(), size0, nBinsActive);fflush(stdout);
    it->second.resize(size0 + nBinsActive);
    offset = size0;
  }
  auto& vec2 = it->second;
  //printf("offset = %lu\n", offset);fflush(stdout);
  for(int i = 0; i < nBinsActive; i++) {
    vec2[offset + i] = vec[i];
  }
}

//
void ReactionFFABM_DISNC_CC::atIteration()
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
    //(_f2abm[ds])[0] = -100.;
    //(_flabm[ds])[0] = -100.;
    //(_f3abm[ds])[0] = -100.;
    _calcf2fldone[ds] = false;
  }

  auto td = _tdDS.begin()->second;
  //td->actualizeWrappers();
  //abm::pdffillgrid();
  _need_pdffillgrid = true;

  // parallel computaion for groups of data points
  //printf("atIteration ReactionTheory::_ncpu = %d _flagComputeAtIteration = %d\n", ReactionTheory::_ncpu, _flagComputeAtIteration);fflush(stdout);
  if (ReactionTheory::_ncpu > 1 || _flagComputeAtIteration) {
    auto td = _tdDS.begin()->second;
    td->actualizeWrappers();
    //printf("ReactionFFABM_DISNC_CC pdffillgrid()\n");
    abm::pdffillgrid();
    _need_pdffillgrid = false;
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
    for (auto ds : _dsIDs) {
      _calcf2fldone[ds] = true;
    }
  }
}

// Place calculations in one function, to optimize calls.
void ReactionFFABM_DISNC_CC::calcF2FL(unsigned dataSetID)
{
  //printf("calcF2FL\n");fflush(stdout);
  //if ((_f2abm[dataSetID][0] < -99.))
  //if ((_f2abm[dataSetID][0] == 0.))
  if(!_calcf2fldone[dataSetID])
  { // compute
    //printf("calcF2FL comp\n");fflush(stdout);
    auto td = _tdDS[dataSetID];
    bool need_update = td->actualizeWrappers();
    if (_need_pdffillgrid || need_update) {
      //printf("ReactionFFABM_DISNC_CC pdffillgrid()\n");
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
    //const size_t np = xp->size();

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
      //for (size_t i = 0; i < np; i++) {
        //auto& point = vec[i];
      for (auto& v : vec) {
        v.calc();
        _f2abm[dataSetID][v.i] = v.f2;
        _flabm[dataSetID][v.i] = v.fl;
        _f3abm[dataSetID][v.i] = v.f3;
      }
    }
    else {
      ForkPool pool(ncpu, _task_distr);
      size_t np = vec.size();
      ForkPool::SharedMemory shm(sizeof(double) * 3 * np);
      double* f2 = shm.data<double>();
      double* fl = f2 + np;
      double* f3 = fl + np;
      pool.parallel_for(vec, [&vec, f2, fl, f3](size_t i) {
        auto& point = vec[i];
        point.calc();
        f2[i] = point.f2;
        fl[i] = point.fl;
        f3[i] = point.f3;
      });
      size_t icomp = 0;
      for (size_t i = 0; i < xp->size(); i++) {
        if(q2[i] > _q2mincomp) {
          _f2abm[dataSetID][i] = f2[icomp];
          _flabm[dataSetID][i] = fl[icomp];
          _f3abm[dataSetID][i] = f3[icomp];
          icomp++;
        }
      }
    }
    _calcf2fldone[dataSetID] = true;
  }
}

// Calculates one data point at (Q2,x) and returns values f2, fl, f3
void ReactionFFABM_DISNC_CC::DataPoint::calc()
{
  static constexpr abm::SFproc ncflag = abm::SFproc::cc;
  //printf("point::calc() x,q2  %f,%f\n", x, q2);fflush(stdout);
  f2 = fl = f3 = 0.;
  double f2l(0), f2b(0), f2c(0), fll(0), flc(0), flb(0), f3l(0), f3b(0), f3c(0);
  if (flav == ReactionBaseDISCC::dataFlav::incl || flav == ReactionBaseDISCC::dataFlav::l) {
    //printf("q2, x, ord, ordHQ, ordFL, msbarmin, charge, sin2thetaW, polar, mz = %f, %f, %d, %d, %d, %d, %f, %f, %f, %f\n", q2, x, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, mz);fflush(stdout);
    f2l = abm::calc_point_strfun(ncflag, abm::SFtype::f2, abm::SFflav::l, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, mz);
    fll = abm::calc_point_strfun(ncflag, abm::SFtype::fl, abm::SFflav::l, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, mz);
    f3l = abm::calc_point_strfun(ncflag, abm::SFtype::f3, abm::SFflav::l, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, mz);
    //printf("f2l, fll, f3l = %f, %f, %f\n", f2l, fll, f3l);fflush(stdout);
  }
  if (flav == ReactionBaseDISCC::dataFlav::incl || flav == ReactionBaseDISCC::dataFlav::c) {
    f2c = abm::calc_point_strfun(ncflag, abm::SFtype::f2, abm::SFflav::c, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, mz);
    flc = abm::calc_point_strfun(ncflag, abm::SFtype::fl, abm::SFflav::c, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, mz);
    f3c = abm::calc_point_strfun(ncflag, abm::SFtype::f3, abm::SFflav::c, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, mz);
  }
  if (flav == ReactionBaseDISCC::dataFlav::incl || flav == ReactionBaseDISCC::dataFlav::b) {
    f2b = abm::calc_point_strfun(ncflag, abm::SFtype::f2, abm::SFflav::b, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, mz);
    flb = abm::calc_point_strfun(ncflag, abm::SFtype::fl, abm::SFflav::b, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, mz);
    f3b = abm::calc_point_strfun(ncflag, abm::SFtype::f3, abm::SFflav::b, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, mz);
  }
  //printf("f2l, f2c, f2b = %f, %f, %f\n", f2l, f2c, f2b);fflush(stdout);
                
  double f3lout_bar = 0.;
  if(nuke && nuke->need_f3bar()) {
    // need F3bar for nuclear corrections and antineutrino
    // we can calculate it now, because HT and TMC (calculated later) do not apply to F3
    int charge_bar = -1 * charge;
    double f2l_bar(0), f2b_bar(0), f2c_bar(0), fll_bar(0), flc_bar(0), flb_bar(0), f3l_bar(0), f3c_bar(0), f3b_bar(0);
    if (flav == ReactionBaseDISCC::dataFlav::incl || flav == ReactionBaseDISCC::dataFlav::l) {
      f3l_bar = abm::calc_point_strfun(ncflag, abm::SFtype::f3, abm::SFflav::l, q2, x, -1, ord, ordHQ, ordHQ, msbarmin, charge_bar, *sin2thetaWPtr, polar, mz);
    }
    if (flav == ReactionBaseDISCC::dataFlav::incl || flav == ReactionBaseDISCC::dataFlav::c) {
      f3c_bar = abm::calc_point_strfun(ncflag, abm::SFtype::f3, abm::SFflav::c, q2, x, -1, ord, ordHQ, ordHQ, msbarmin, charge_bar, *sin2thetaWPtr, polar, mz);
    }
    if (flav == ReactionBaseDISCC::dataFlav::incl || flav == ReactionBaseDISCC::dataFlav::b) {
      f3b_bar = abm::calc_point_strfun(ncflag, abm::SFtype::f3, abm::SFflav::b, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge_bar, *sin2thetaWPtr, polar, mz);
    }
    f3lout_bar = x * ReactionFFABM_DISNC_CC::combine_flavours(flav, f3l_bar, f3c_bar, f3b_bar);
  }
  if (tmc) {
    static constexpr abm::SFproc ncflag = abm::SFproc::nc;
    tmc->apply(f2l, fll, f3l, f2c, flc, f3c, f2b, flb, f3b, q2, x, ncflag, ord, ordHQ, ordFL, msbarmin, charge, polar, 1.-*sin2thetaWPtr, mz);
  }
  if(ht) {
    // HT is applied only to F2 and FL light flavour part
    ht->apply(q2, x, f2l, fll, f3l);
  }
  //printf("f2l, f2c, f2b = %f, %f, %f\n", f2l, f2c, f2b);fflush(stdout);
  f2 = ReactionFFABM_DISNC_CC::combine_flavours(flav, f2l, f2c, f2b);
  fl = ReactionFFABM_DISNC_CC::combine_flavours(flav, fll, flc, flb);
  f3 = x * combine_flavours(flav, f3l, f3c, f3b);
  // apply nuclear corrections to the sum of light+c+b because corrections for charm and non-charm (kint=4,5) are not implemented
  if (nuke) {
    nuke->apply(q2, x, f2, fl, f3);
  }
  //printf("f2, fl, f3 = %f, %f, %f %p\n", f2, fl, f3, this);fflush(stdout);
}

double ReactionFFABM_DISNC_CC::combine_flavours(const ReactionBaseDISCC::dataFlav flav, const double f, const double fc, const double fb)
{
  switch (flav)
  {
    case ReactionBaseDISCC::dataFlav::incl:
      return f + fc + fb;
    case ReactionBaseDISCC::dataFlav::c:
      return fc;
    case ReactionBaseDISCC::dataFlav::b:
      return fb;
    default:
      hf_errlog(28022501, "F: Unsupported flavour");
      return 0.; // avoid warning
  }
}

valarray<double> ReactionFFABM_DISNC_CC::F2(TermData *td)
{
  calcF2FL(td->id);
  return _f2abm[td->id];
}

valarray<double> ReactionFFABM_DISNC_CC::FL(TermData *td)
{
  calcF2FL(td->id);
  return _flabm[td->id];
}

valarray<double> ReactionFFABM_DISNC_CC::xF3(TermData *td)
{
  calcF2FL(td->id);
  return _f3abm[td->id];
}
