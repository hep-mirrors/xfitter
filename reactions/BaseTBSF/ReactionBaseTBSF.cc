#include "ReactionBaseTBSF.h"
#include "ReactionBaseDISNC.h"
#include "ReactionBaseDISCC.h"
/*#include "DIS_HT.h"
#include "DIS_TMC.h"
#include "DIS_NUKE.h"
#include "TermData.h"
#include "ForkPool.h"
#include "xfitter_cpp_base.h"
#include "xfitter_steer.h"

extern "C" {
  double numufcalflux_(const double& e); // NOMAD E(nu) flux
  double sd2_(double* acc, double (*f)(double*), void (*r)(int, double*, double*, double*)); // integrator
}

template <class BaseDIS, abm::SFproc Proc>
ReactionBaseFFABM<BaseDIS, Proc>::~ReactionBaseFFABM() {
  for (auto ht : _ht) delete ht.second;
  for (auto tmc : _tmc) delete tmc.second;
  for (auto nk : _nuke) delete nk.second;
  //for (const auto& it : _grouped_shared_memory) delete it.second;
}

template <class BaseDIS, abm::SFproc Proc>
void ReactionBaseFFABM<BaseDIS, Proc>::initTerm(TermData *td) {
  BaseDIS::initTerm(td);
  unsigned termID = td->id;

  // scales mu^2 = scalea1 * Q^2 + scaleb1 * 4 * m_h^2 (default scalea1 = scaleb1 = 1.0)
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

  // CKM matrix
  _ckm.resize(9);
  _ckm[0] = td->getParamD("Vud");
  _ckm[1] = td->getParamD("Vus");
  _ckm[2] = td->getParamD("Vub");
  _ckm[3] = td->getParamD("Vcd");
  _ckm[4] = td->getParamD("Vcs");
  _ckm[5] = td->getParamD("Vcb");
  _ckm[6] = td->getParamD("Vtd");
  _ckm[7] = td->getParamD("Vts");
  _ckm[8] = td->getParamD("Vtb");

  // read higher twist parameters
  _ht[termID] = nullptr;
  if (td->hasParam("ht") && td->getParamI("ht") != 0) {
    _ht[termID] = new DIS_HT(td);
  }
  // read target mass correction parameters
  _tmc[termID] = nullptr;
  if (td->hasParam("tmc") && td->getParamS("tmc") != "") {
    _tmc[termID] = new DIS_TMC(td);
  }
  // read nuclear correction parameters
  _nuke[termID] = nullptr;
  if (td->hasParam("nuke_ftyp") && td->getParamI("nuke_ftyp") != 0) {
    _nuke[termID] = new DIS_NUKE(td);
  }

  const std::valarray<double>* q2_ptr = nullptr;
  const std::valarray<double>* x_ptr = nullptr;
  const std::valarray<double>* onedimvar_ptr = nullptr;
  const std::valarray<double>* q2min_ptr = nullptr;
  const std::valarray<double>* q2max_ptr = nullptr;
  const std::valarray<double>* ymin_ptr = nullptr;
  const std::valarray<double>* ymax_ptr = nullptr;
  double energy;
  int nBins;
  auto find_binning = [this,td,&nBins,&q2_ptr,&x_ptr,&onedimvar_ptr,&q2min_ptr,&q2max_ptr,&ymin_ptr,&ymax_ptr,&energy]() {
    if (td->hasParam("binning")) {
      auto binning = td->getParamS("binning");
      if(binning == "E") {
        onedimvar_ptr = td->getBinColumnOrNull("E");
        nBins = onedimvar_ptr->size();
        return Binning::point_at_e;
      }
      else if(binning == "x") {
        onedimvar_ptr = td->getBinColumnOrNull("x");
        nBins = onedimvar_ptr->size();
        return Binning::point_at_x;
      }
      else if(binning == "sqrtshat") {
        onedimvar_ptr = td->getBinColumnOrNull("sqrtshat");
        nBins = onedimvar_ptr->size();
        return Binning::point_at_sqrtshat;
      }
      else if(binning == "bin_q2y") {
        q2min_ptr = td->getBinColumnOrNull("q2min");
        q2max_ptr = td->getBinColumnOrNull("q2max");
        ymin_ptr = td->getBinColumnOrNull("ymin");
        ymax_ptr = td->getBinColumnOrNull("ymax");
        energy = *td->getParamD("energy");
        nBins = q2min_ptr->size();
        return Binning::bin_q2y;
      }
      else {
        hf_errlog(2026062202, "F: Unsupported data point binning " + binning);
      }
    }
    else {
      //q2_ptr = td->getBinColumnOrNull("Q2");
      //x_ptr = td->getBinColumnOrNull("x");
      q2_ptr = GetBinValues(td, "Q2");
      x_ptr = GetBinValues(td, "x");
      if(q2_ptr && x_ptr) {
        nBins = q2_ptr->size();
        return Binning::point_at_q2x;
      }
      onedimvar_ptr = td->getBinColumnOrNull("E");
      if(onedimvar_ptr) {
        nBins = onedimvar_ptr->size();
        return Binning::point_at_e;
      }
      onedimvar_ptr = td->getBinColumnOrNull("x");
      if(onedimvar_ptr) {
        nBins = onedimvar_ptr->size();
        return Binning::point_at_x;
      }
      onedimvar_ptr = td->getBinColumnOrNull("sqrtshat");
      if(onedimvar_ptr) {
        nBins = onedimvar_ptr->size();
        return Binning::point_at_sqrtshat;
      }
      q2min_ptr = td->getBinColumnOrNull("q2min");
      q2max_ptr = td->getBinColumnOrNull("q2max");
      ymin_ptr = td->getBinColumnOrNull("ymin");
      ymax_ptr = td->getBinColumnOrNull("ymax");
      energy = *td->getParamD("energy");
      if(q2min_ptr && q2max_ptr && ymax_ptr) {
        nBins = q2min_ptr->size();
        return Binning::bin_q2y;
      }
    }
    hf_errlog(26062301, "F: cannot determine binning");
    return Binning::bin_q2y; // avoid warning
  };
    Binning datapoint_binning = find_binning();
  const int scalesemilepbr = td->getParamI("scalesemilepbr");
  const double* br0 = td->hasParam("br_cmu_0") ? td->getParamD("br_cmu_0") : nullptr;
  const double* br1 = td->hasParam("br_cmu_1") ? td->getParamD("br_cmu_1") : nullptr;
  const double mnucl = td->hasParam("mpr") ? *td->getParamD("mpr") : 0.0;
  const double mw = *td->getParamD("Mw");
  //const double mnucl = (*td->getParamD("mpr") + *td->getParamD("mnt"))/2.;
  const auto& integrator_str = td->getParamS("integrator");
  Integrator integrator;
  if(integrator_str == "sd2") {
    integrator = Integrator::sd2;
  }
  else {
    hf_errlog(26061101, "F: unknown integrator " + integrator_str);
  }
  const double integrator_epsrel = *td->getParamD("integrator_epsrel");
  const int integrator_verbose = td->getParamI("integrator_verbose");

  printf("---------------------------------------------\n");
  printf("INFO from %s: ", BaseDIS::getReactionName().c_str());
  printf("running mass = %c", msbarmin ? 'T' : 'F');
  printf(", order = %d", order);
  printf(", order HQ = %d", orderHQ);
  printf(", O(alpha_S) F_L - O(alpha_S) F2 = %d", ordfl);
  printf(", factorisation scale for heavy quarks = sqrt(%f * Q^2 + %f * 4m_q^2", hqscale1in, hqscale2in);
  printf(", binning = %d\n", int(datapoint_binning));
  printf("---------------------------------------------\n");

  _f2abm[termID].resize(nBins);
  _flabm[termID].resize(nBins);
  _f3abm[termID].resize(nBins);

  const double* mzPtr = td->getParamD("Mz");
  const double* sin2thwPtr = td->getParamD("sin2thW");

  // groups for parallel
  const auto& group = td->getParamS("group_parallel");
  printf("group = %s\n", group.c_str());
  size_t offset = 0;
  size_t nBinsActive = nBins;
  if(datapoint_binning == Binning::point_at_q2x) {
    nBinsActive = std::count_if(std::begin(*q2_ptr), std::end(*q2_ptr), [this](double q2) { return q2 >= this->_q2mincomp; });
  }
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
    if(datapoint_binning == Binning::point_at_q2x) {
      if((*q2_ptr)[i] <= _q2mincomp) {
        _f2abm[i] = 0.0;
        _flabm[i] = 0.0;
        _f3abm[i] = 0.0;
        continue;
      }
    }
    auto& point = vec2[offset + iActive];
    point.binning = datapoint_binning;
    if(datapoint_binning == Binning::point_at_q2x) {
      point.q2 = (*q2_ptr)[i];
      point.x = (*x_ptr)[i];
    }
    else if(datapoint_binning == Binning::point_at_e || datapoint_binning == Binning::point_at_x || datapoint_binning == Binning::point_at_sqrtshat) {
      point.onedimvar = (*onedimvar_ptr)[i];
    }
    else if(datapoint_binning == Binning::bin_q2y) {
      point.q2min = (*q2min_ptr)[i];
      point.q2max = (*q2max_ptr)[i];
      if(ymin_ptr) {
        point.ymin = (*ymin_ptr)[i];
      }
      else {
        point.ymin = 0.0;
      }
      point.ymax = (*ymax_ptr)[i];
      point.energy = energy;
    }
    point.mnucl = mnucl;
    point.mw = mw;
    point.scalesemilepbr = scalesemilepbr;
    point.br0 = br0;
    point.br1 = br1;
    point.integrator = integrator;
    point.integrator_epsrel = integrator_epsrel;
    point.integrator_verbose = integrator_verbose;
    point.td = td;
    point.i = i;
    point.is_beam_nu = BaseDIS::IsBeamNu(td->id);
    point.flav = dataFlav(BaseDIS::GetDataFlav(td->id));
    point.ord = order;
    point.ordHQ = orderHQ;
    point.ordFL = ordfl;
    point.msbarmin = msbarmin;
    point.charge = BaseDIS::GetCharge(td->id);
    point.polar = BaseDIS::GetPolarisation(td->id);
    point.sin2thetaWPtr = sin2thwPtr;
    point.mz = mzPtr;
    point.ht = _ht[td->id];
    point.tmc = _tmc[td->id];
    point.nuke = _nuke[td->id];
    point.f2 = 0.0;
    point.fl = 0.0;
    point.f3 = 0.0;
    iActive++;
  }
}

template <class BaseDIS, abm::SFproc Proc>
void ReactionBaseFFABM<BaseDIS, Proc>::atIteration() {
  BaseDIS::atIteration();
  abm::set_hq_masses(*_mcPtr, *_mbPtr);
  if(Proc == abm::SFproc::cc) {
    abm::update_ckm_matrix(_ckm);
  }
  for (auto ht : _ht) {
    if (ht.second) {
      ht.second->update();
    }
  }
  _grouped_data_points.begin()->second[0].td->actualizeWrappers();
  abm::pdffillgrid();

  // parallel computaion for groups of data points
  //printf("atIteration _ncpu = %d\n", BaseDIS::_ncpu);
  int ncpu =  xfitter::xf_ncpu(BaseDIS::_ncpu);
  if(ncpu == 1) {
    for(auto& it : _grouped_data_points) {
      for(auto& point : it.second) {
        point.calc();
        _f2abm[point.td->id][point.i] = point.f2;
        _flabm[point.td->id][point.i] = point.fl;
        _f3abm[point.td->id][point.i] = point.f3;
      }
    }
  }
  else {
    for(auto& it : _grouped_data_points) {
      auto& vec = it.second;
      bool need_pdffillgrid = vec[0].td->actualizeWrappers();
      printf("checking pdffillgrid() group = %s, datasetID = %d\n", it.first.c_str(), vec[0].td->id);
      if(need_pdffillgrid) {
        printf("extra pdffillgrid() group = %s, datasetID = %d\n", it.first.c_str(), vec[0].td->id);
        abm::pdffillgrid();
      }
      ForkPool pool(ncpu, BaseDIS::_task_distr);
      int np = vec.size();
      ForkPool::SharedMemory shm(sizeof(double) * 3 * np + sizeof(int) * np * 2);
      //if(_grouped_shared_memory.find(it.first) == _grouped_shared_memory.end()) {
      //  _grouped_shared_memory[it.first] = new ForkPool::SharedMemory(sizeof(double) * 3 * np + sizeof(int) * np * 2);
      //}
      //ForkPool::SharedMemory& shm = *_grouped_shared_memory[it.first];
      double* f2 = shm.data<double>();
      double* fl = f2 + np;
      double* f3 = fl + np;
      int* datasetID = (int*)(f3 + np);
      int* i_orig = datasetID + np;
      pool.parallel_for(np, [&vec, f2, fl, f3, datasetID, i_orig](size_t i) {
        auto& point = vec[i];
        point.calc();
        f2[i] = point.f2;
        fl[i] = point.fl;
        f3[i] = point.f3;
        datasetID[i] = point.td->id;
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

template <class BaseDIS, abm::SFproc Proc>
double ReactionBaseFFABM<BaseDIS, Proc>::DataPoint::integrand_sd2(double x[]) {
  const int ndim = 2;
  const int ncomp = 1;
  double val = 0.0;
  integrand(&ndim, x, &ncomp, &val, &_integration_params_static);
  return val;
}

template <class BaseDIS, abm::SFproc Proc>
int ReactionBaseFFABM<BaseDIS, Proc>::DataPoint::integrand(const int* ndim, const double* inp, const int *ncomp, double* val, void *pars_voidptr) {
  //(void)ndim; // Unused parameter
  static constexpr double xmax = 0.99;
  static constexpr double q2min = 1.;
  integration_params& pars = *(integration_params*)pars_voidptr;
  pars.ncalls++;
  double x(-1.), q2(-1.), y(-1.), e(-1.), s(-1.), factor(-1.);
  switch(pars.point->binning) {
    case Binning::point_at_e: {
      e = pars.point->onedimvar;
      s = 2 * pars.point->mnucl * e;
      double xmin1 = q2min / s;
      x = xmin1 + (xmax - xmin1) * inp[0];
      double ymin = q2min / s / x;
      y = ymin + (1 - ymin) * inp[1];
      q2 = s * x * y;
      factor = (xmax - xmin1) * (1 - ymin);
      break;
    }
    case Binning::point_at_x: {
      x = pars.point->onedimvar;
      double alim = q2min / x / 2. / pars.point->mnucl;
      double emin1 = std::max(std::min(alim, pars.emax), pars.emin);
      e = emin1 + (pars.emax - emin1) * inp[0];
      double ymin1 = alim / e;
      y = ymin1 + (1 - ymin1) * inp[1];
      s = 2 * pars.point->mnucl * e;
      q2 = s * x * y;
      factor = pars.eflux(e);
      factor *= (pars.emax - emin1) * (1 - ymin1);
      factor *= e;
      break;
    }
    case Binning::point_at_sqrtshat: {
      double sqrtshat = pars.point->onedimvar;
      double shat = sqrtshat*sqrtshat;
      double emin1 = std::min(std::max(pars.emin, (q2min + shat) / 2. / pars.point->mnucl), pars.emax);
      e = emin1 + (pars.emax - emin1) * inp[0];
      double spmax = 2 * pars.point->mnucl * e;
      double xmin1 = 1. / (1 + shat / q2min);
      double xmax1 = 1 - shat / spmax;
      x = xmin1 + (xmax1 - xmin1) * inp[1];
      q2 = shat * x / (1 - x);
      y = q2 / (2 * e * pars.point->mnucl * x);
      factor = pars.eflux(e);
      factor *= 1. / (2 * pars.point->mnucl * (1 - x));
      factor *= (xmax1 - xmin1) * (pars.emax - emin1);
      break;
    }
    case Binning::bin_q2y: {
      s = pars.point->energy * pars.point->energy;
      q2 = pars.point->q2min + (pars.point->q2max - pars.point->q2min) * inp[0];
      double ymin1 = std::max(pars.point->q2min / s, pars.point->ymin);
      y = ymin1 + (1 - ymin1) * inp[1];
      x = q2 / s / y;
      //printf("q2,x,y = %f,%f,%f\n", q2,x,y);
      factor = (pars.point->q2max - pars.point->q2min) * (1 - ymin1);
      break;
    }
    default:
      hf_errlog(26062302, "F: Unsupported integration with binning " + std::to_string(int(pars.point->binning)));
      return 0.; // avoid warning
  }
  if (x <= 0. || x >= 1.) {
    val[0] = 0.;
    return 0;
  }
  if (y <= 0. || y >= 1.) {
    val[0] = 0.;
    return 0;
  }
  pars.point->q2 = q2;
  pars.point->x = x;
  pars.point->calc_at_q2x();
  double yplus = 1.0 + (1.0 - y) * (1.0 - y) - 2.0 * std::pow(pars.point->mnucl * x * y, 2.0) / q2;
  double yminus = 1.0 - (1.0 - y) * (1.0 - y);
  auto charge_mod = pars.point->is_beam_nu ? -1 * pars.point->charge : pars.point->charge; // swap charge for nu beam, see InitTerm()
  val[0] = 0.5 * (1 + charge_mod * pars.point->polar) * (yplus * pars.point->f2 - charge_mod * yminus * pars.point->f3 - y * y * pars.point->fl);
  val[0] *= factor;
  if (pars.point->scalesemilepbr) {
    if(!pars.point->br0 || !pars.point->br1) {
      hf_errlog(2026070201, "F: Parameters br_cmu_0 or br_cmu_1 were not provided");
    }
    double br = *pars.point->br0 / (1 + *pars.point->br1 / e);
    val[0] *= br;
  }
  if (pars.point->flav == dataFlav::incl || pars.point->flav == dataFlav::l) {
    val[0] *= std::pow((pars.point->mw * pars.point->mw / (q2 + pars.point->mw * pars.point->mw)), 2.);
  }
  //if (pars.nomad_scaleq2mw2) {
  //  if (pars.rd->_dataFlav == BaseDISCC::dataFlav::incl || pars.rd->_dataFlav == BaseDISCC::dataFlav::l) {
  //    val[0] *= std::pow((MW * MW / (q2 + MW * MW)), 2.);
  //  }
  //}
  if (val[0] != val[0] || 1. / val[0] == 0.) {
    val[0] = 0.;
  }
  return 0;
}

template <class BaseDIS, abm::SFproc Proc>
void ReactionBaseFFABM<BaseDIS, Proc>::DataPoint::integrand_sd2_region(int ll, double* xx, double* aa, double* bb)
{
  (void)ll;
  (void)xx;
  const double del = 1e-8;
  //const double del = 0.0;
  aa[0] = del;
  bb[0] = 1.0 - del;
  aa[1] = del;
  bb[1] = 1.0 - del;
}

template <class BaseDIS, abm::SFproc Proc>
void ReactionBaseFFABM<BaseDIS, Proc>::DataPoint::calc() {
  if(binning == Binning::point_at_q2x) {
    calc_at_q2x();
  }
  else if(binning == Binning::point_at_e || binning == Binning::point_at_x || binning == Binning::point_at_sqrtshat || binning == Binning::bin_q2y) {
    calc_2d_integral();
  }
  else {
    hf_errlog(2026062203, "F: Unsupported data point binning " + std::to_string(int(binning)));
  }
}

template <class BaseDIS, abm::SFproc Proc>
void ReactionBaseFFABM<BaseDIS, Proc>::DataPoint::calc_2d_integral() {
  _integration_params_static.point = this;
  _integration_params_static.eflux = numufcalflux_;
  _integration_params_static.emin = 6.;
  _integration_params_static.emax = 300.;
  _integration_params_static.ncalls = 0;
  switch(integrator) {
    case Integrator::sd2: {
      double acc = 0.;
      double s0=sd2_(&acc, integrand_sd2, integrand_sd2_region);
      acc = s0 * integrator_epsrel;
      double s1=sd2_(&acc, integrand_sd2, integrand_sd2_region);
      if(integrator_verbose) {
        printf("binning,onedimvar = %d,%e s0,s1 = %e,%e\n", int(_integration_params_static.point->binning), _integration_params_static.point->onedimvar, s0, s1);
      }
      f2 = s1;
      break;
    }
  }
  if(integrator_verbose) {
    printf("ncalls = %ld\n", _integration_params_static.ncalls);
  }
}

template <class BaseDIS, abm::SFproc Proc>
void ReactionBaseFFABM<BaseDIS, Proc>::DataPoint::calc_at_q2x() {
  auto combine_flavours = [](const dataFlav flav, const double l, const double c, const double b) {
    switch (flav)
    {
      case dataFlav::incl:
        return l + c + b;
      case dataFlav::l:
        return l;
      case dataFlav::c:
        return c;
      case dataFlav::b:
        return b;
      default:
        hf_errlog(28022501, "F: Unsupported flavour");
        return 0.; // avoid warning
    }
  };

  bool need_pdffillgrid = td->actualizeWrappers();
  if(need_pdffillgrid) {
    printf("extra pdffillgrid() in calc_at_q2x() by datasetID = %d point = %d (proc_NCCC = %d)\n", td->id, i, int(Proc));
    //fflush(stdout);
    abm::pdffillgrid();
  }
  printf("q2,x = %f,%f\n", q2, x);fflush(stdout);
  f2 = fl = f3 = 0.;
  double f2l(0), f2b(0), f2c(0), fll(0), flc(0), flb(0), f3l(0), f3b(0), f3c(0);
  if (flav == dataFlav::incl || flav == dataFlav::l) {
    f2l = abm::calc_point_strfun(Proc, abm::SFtype::f2, abm::SFflav::l, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, *mz);
    fll = abm::calc_point_strfun(Proc, abm::SFtype::fl, abm::SFflav::l, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, *mz);
    f3l = abm::calc_point_strfun(Proc, abm::SFtype::f3, abm::SFflav::l, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, *mz);
  }
  if (flav == dataFlav::incl || flav == dataFlav::c) {
    f2c = abm::calc_point_strfun(Proc, abm::SFtype::f2, abm::SFflav::c, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, *mz);
    flc = abm::calc_point_strfun(Proc, abm::SFtype::fl, abm::SFflav::c, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, *mz);
    f3c = abm::calc_point_strfun(Proc, abm::SFtype::f3, abm::SFflav::c, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, *mz);
  }
  if (flav == dataFlav::incl || flav == dataFlav::b) {
    f2b = abm::calc_point_strfun(Proc, abm::SFtype::f2, abm::SFflav::b, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, *mz);
    flb = abm::calc_point_strfun(Proc, abm::SFtype::fl, abm::SFflav::b, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, *mz);
    f3b = abm::calc_point_strfun(Proc, abm::SFtype::f3, abm::SFflav::b, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge, *sin2thetaWPtr, polar, *mz);
  }
                
  double f3lout_bar = 0.;
  if(nuke && nuke->need_f3bar()) {
    // need F3bar for nuclear corrections and antineutrino
    // we can calculate it now, because HT and TMC (calculated later) do not apply to F3
    int charge_bar = -1 * charge;
    double f2l_bar(0), f2b_bar(0), f2c_bar(0), fll_bar(0), flc_bar(0), flb_bar(0), f3l_bar(0), f3c_bar(0), f3b_bar(0);
    if (flav == dataFlav::incl || flav == dataFlav::l) {
      f3l_bar = abm::calc_point_strfun(Proc, abm::SFtype::f3, abm::SFflav::l, q2, x, -1, ord, ordHQ, ordHQ, msbarmin, charge_bar, *sin2thetaWPtr, polar, *mz);
    }
    if (flav == dataFlav::incl || flav == dataFlav::c) {
      f3c_bar = abm::calc_point_strfun(Proc, abm::SFtype::f3, abm::SFflav::c, q2, x, -1, ord, ordHQ, ordHQ, msbarmin, charge_bar, *sin2thetaWPtr, polar, *mz);
    }
    if (flav == dataFlav::incl || flav == dataFlav::b) {
      f3b_bar = abm::calc_point_strfun(Proc, abm::SFtype::f3, abm::SFflav::b, q2, x, -1, ord, ordHQ, ordFL, msbarmin, charge_bar, *sin2thetaWPtr, polar, *mz);
    }
    f3lout_bar = x * combine_flavours(flav, f3l_bar, f3c_bar, f3b_bar);
  }
  if (tmc) {
    tmc->apply(f2l, fll, f3l, f2c, flc, f3c, f2b, flb, f3b, q2, x, Proc, ord, ordHQ, ordFL, msbarmin, charge, polar, 1.-*sin2thetaWPtr, *mz);
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
    nuke->apply(q2, x, f2, fl, f3, &f3lout_bar);
  }
}

template <class BaseDIS, abm::SFproc Proc>
typename ReactionBaseFFABM<BaseDIS, Proc>::integration_params ReactionBaseFFABM<BaseDIS, Proc>::DataPoint::_integration_params_static;
*/

template class ReactionBaseTBSF<ReactionBaseDISNC, abm::SFproc::nc>;
template class ReactionBaseTBSF<ReactionBaseDISCC, abm::SFproc::cc>;
