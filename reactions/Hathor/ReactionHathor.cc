
/*
   @file ReactionHathor.cc
   @date 2017-08-07
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-08-07
*/

#include "ReactionHathor.h"
#include "HathorPdfxFitter.h"
#include "Hathor.h"
#include "HathorGenericIntegrator.h"
#include "cstring"
#include "xfitter_cpp.h"
#include "xfitter_steer.h"
#include <unistd.h>

// the class factories
extern "C" ReactionHathor* create()
{
  return new ReactionHathor();
}

extern "C"
{
int rlxd_size(void);
void rlxd_get(int state[]);
void rlxd_reset(int state[]);
void rlxd_init(int level,int seed);
}

ReactionHathor::ReactionHathor()
{
  _pdf = NULL;
  _rndStore = NULL;
  _hathor = NULL;
  //_mtop = -1.0;
  //_mr = -1.0;
  //_mf = -1.0;
}

ReactionHathor::~ReactionHathor()
{
  if(_rndStore)
    delete[] _rndStore;

  // do NOT delete Hathor instances here, because:
  // (1) Hathor classes do not have vitual destructors which produces a warning
  // (2) Hathor classes do not have destructors at all

  //if(_pdf)
  //  delete _pdf;
  //for(auto item : _hathorArray)
  //  if(item.second)
  //    delete item.second;
}

void ReactionHathor::atStart()
{
  ReactionTheory::atStart();

  // PDFs for Hathor
  _pdf = new HathorPdfxFitter();

    // random number generator
  rlxd_init(1, 1);
  int nRnd = rlxd_size();
  //std::cout << " Size of random number array = " << nRnd << "\n";
  _rndStore = new int[nRnd];
  rlxd_get(_rndStore);

  // instantiate one Hathor instance for all terms
  //_hathor = new Hathor(*_pdf);
  _hathor = new HathorGenericIntegrator<HathorLikeSgTop>(*_pdf);
}

void ReactionHathor::initTerm(TermData *td)
{
  ReactionTheory::initTerm(td);
  int dataSetID = td->id;

  // check if dataset with provided ID already exists
  if(_mtopPerInstance.find(dataSetID) != _mtopPerInstance.end() && !chi2scan_.scan)
  {
    char str[256];
    sprintf(str, "F: dataset with id = %d already exists", dataSetID);
    hf_errlog_(17080701, str, strlen(str));
  }

  _integrator = td->getParamS("integrator");
  _mtopPerInstance[dataSetID] = td->getParamD("mtp");
  _mfPerInstance[dataSetID] = *td->getParamD("muR");
  _mrPerInstance[dataSetID] = *td->getParamD("muF");
  _orderPerInstance[dataSetID] = td->getParamS("Order");
  _NfPerInstance[dataSetID] = td->hasParam("NFlavour") ? td->getParamI("NFlavour") : 5;
  _MsbarPerInstance[dataSetID] = td->getParamI("MS_MASS");
  _precisionLevelPerInstance[dataSetID] = td->getParamI("precisionLevel");
  _sqrtSPerInstance[dataSetID] = *td->getParamD("SqrtS");
  _ppbarPerInstance[dataSetID] = td->getParamI("ppbar");


  // avoid calculating the same cross section many times
  std::string calc_name = "mtp" + std::to_string((unsigned long long)(void**)_mtopPerInstance[dataSetID]);
  calc_name += std::string("_") + "muF" + std::to_string(_mfPerInstance[dataSetID]);
  calc_name += std::string("_") + "muR" + std::to_string(_mrPerInstance[dataSetID]);
  calc_name += std::string("_") + "Order" + _orderPerInstance[dataSetID];
  calc_name += std::string("_") + "NFlavour" + std::to_string(_NfPerInstance[dataSetID]);
  calc_name += std::string("_") + "MS_MASS" + std::to_string(_MsbarPerInstance[dataSetID]);
  calc_name += std::string("_") + "precisionLevel" + std::to_string(_precisionLevelPerInstance[dataSetID]);
  calc_name += std::string("_") + "SqrtS" + std::to_string(_sqrtSPerInstance[dataSetID]);
  calc_name += std::string("_") + "ppbar" + std::to_string(_ppbarPerInstance[dataSetID]);
  if(td->hasParam("evolution")) {
    calc_name += "_evolution_" + td->getParamS("evolution");
  }
  if(td->hasParam("evolution1")) {
    calc_name += "_evolution1_" + td->getParamS("evolution1");
  }
  if(td->hasParam("evolution2")) {
    calc_name += "_evolution2_" + td->getParamS("evolution2");
  }
  if (_convolved.find(calc_name) == _convolved.end()) {
    //printf("adding calc_name = %s\n", calc_name.c_str());
    _convolved[calc_name] = std::make_pair<double, std::vector<TermData*> >(0., {});
    _convolved_vector_of_keys.push_back(calc_name);
  }
  _convolved[calc_name].second.push_back(td);
}

// Main function to compute results at an iteration
void ReactionHathor::compute(TermData *td, valarray<double> &val, map<string, valarray<double> > &err)
{
  for(const auto& it : _convolved) {
    for(const auto& p : it.second.second) {
      if(td == p) {
        val = it.second.first;
        return;
      }
    }
  }
  hf_errlog(2026062401, "F: HATHOR no computed cross section for datasetID = " + std::to_string(td->id));
}

void ReactionHathor::atIteration() {
  ReactionTheory::atIteration();

  auto calc_one = [&](int i, double& xsec) {
    TermData* td = _convolved[_convolved_vector_of_keys[i]].second[0];
    td->actualizeWrappers();
    int dataSetID = td->id;
    _pdf->IsValid = true;
    rlxd_reset(_rndStore);

    //Suppress Hathor output
    // SZ 2023.11.11 freopen breaks pipe redicrection, see https://c-faq.com/stdio/undofreopen.html
    // Solution from https://stackoverflow.com/questions/1908687/how-to-redirect-the-output-back-to-the-screen-after-freopenout-txt-a-stdo
    int o;
    if (!steering_.ldebug && _init) {
      o = dup(fileno(stdout));
      if (!freopen("/dev/null", "a", stdout)) {
        hf_errlog(2026062501, "F: Failed to redirect stdout to /dev/null");
      }
    }

    // collision type
    if(_ppbarPerInstance[dataSetID])
      _hathor->setColliderType(Hathor::PPBAR);
    else
      _hathor->setColliderType(Hathor::PP);

    // read centre-of-mass energy
    _hathor->setSqrtShad(_sqrtSPerInstance[dataSetID]);

    // set conversion factor
    static constexpr double convFac_in = 0.389379323e9;  //HATHOR default value
    _hathor->sethc2(convFac_in);

    // set perturbative order and mass scheme
    int scheme = Hathor::LO;
    if(_orderPerInstance[dataSetID] == "NLO")
      scheme = scheme | Hathor::NLO;
    else if(_orderPerInstance[dataSetID] == "NNLO")
      scheme = scheme | Hathor::NLO | Hathor::NNLO;
    if(_MsbarPerInstance[dataSetID])
      scheme = scheme | Hathor::MS_MASS;
    _hathor->setScheme(scheme);

    // number of massless flavours
    _hathor->setNf(_NfPerInstance[dataSetID]);

    // set precision level
    int precisionLevel = std::pow(10, 2 + _precisionLevelPerInstance[dataSetID]);
    // check that this setting is allowed
    // see in AbstractHathor.h:
    //   enum ACCURACY { LOW=1000, MEDIUM=10000, HIGH=100000 };
    // and
    // precisionLevel = 1 -> Hathor::LOW
    // precisionLevel = 2 -> Hathor::MEDIUM
    // precisionLevel = 3 -> Hathor::HIGH
    if(_integrator == "vegas" && precisionLevel !=  Hathor::LOW && precisionLevel !=  Hathor::MEDIUM && precisionLevel !=  Hathor::HIGH) {
      hf_errlog(17081102, "F: provided precision level = " + std::to_string(precisionLevel) + " not supported by Hathor");
    }
    _hathor->setPrecision(precisionLevel);

    // quark mass and scales
    double mt = *_mtopPerInstance[dataSetID];
    double mr = _mrPerInstance[dataSetID] * mt;
    double mf = _mfPerInstance[dataSetID] * mt;
    if (mr<0.) {
      mr = mr*-1./mt;
    }
    if (mf<0.) {
      mf = mf*-1./mt;
    }

    if (steering_.ldebug || !_init)
      {
        std::cout << " Hathor will use for this instance (" + std::to_string(dataSetID) + "):";
        std::cout << " mtop = " << mt << "[GeV] ";
        std::cout << " renorm. scale = " << mr << "[GeV] ";
        std::cout << " factor. scale = " << mf << "[GeV] ";
        std::cout << " SqrtS = " << _sqrtSPerInstance[dataSetID];
        std::cout << " scheme: " << scheme;
        std::cout << " precisionLevel: " << precisionLevel << std::endl;

        // done
        _hathor->PrintOptions();
      }

    _hathor->getXsection(mt, mr, mf, td->getParamS("integrator"));
    double dum = 0.0;
    _hathor->getResult(0, xsec, dum);

    //Resume standard output
    if (!steering_.ldebug && _init)
    {
      dup2(o,fileno(stdout));
      close(o);
    }

    //printf("getXsection ncalls = %ld\n", _ncalls);
    //printf("mt,mr,mf,xsec,err: %f %f %f %f %f [%.3f%%]\n", mt, mr, mf, xsec, dum, dum/xsec*100.);
  };

  int ncpu =  xfitter::xf_ncpu(_ncpu);
  //int ncpu = _ncpu;
  size_t n = _convolved.size();
  if (ncpu == 1) {
    for (int i = 0; i < n; i++) {
      calc_one(i, _convolved[_convolved_vector_of_keys[i]].first);
    }
  }
  else {
    ForkPool pool(ncpu, _task_distr);
    ForkPool::SharedMemory shm(sizeof(double) * n);
    double* val = shm.data<double>();
    pool.parallel_for(n, [&](size_t i) {
      calc_one(i, val[i]);
    });
    for (size_t i = 0; i < n; i++) {
      _convolved[_convolved_vector_of_keys[i]].first = val[i];
    }
  }
  _init = true;
}
