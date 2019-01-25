
/*
   @file ReactionHathor.cc
   @date 2017-08-07
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-08-07
*/

#include "ReactionHathor.h"
#include "HathorPdfxFitter.h"
#include "Hathor.h"
#include "cstring"

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

int ReactionHathor::initAtStart(const string &s)
{
  // PDFs for Hathor
  _pdf = new HathorPdfxFitter(this);

  // random number generator
  rlxd_init(1, 1);
  int nRnd = rlxd_size();
  //std::cout << " Size of random number array = " << nRnd << "\n";
  _rndStore = new int[nRnd];
  rlxd_get(_rndStore);

  // top quark mass
  //std::string mtopName = "mtp";// shouldn't we distinguish somehow between pole and running masses?
  //if(!checkParam(mtopName))
  //{
  //  std::string str = "F: no top quark mass (\"" + mtopName + "\" parameter) for Hathor";
  //  hf_errlog_(17081101, str.c_str(), strlen(str.c_str()));
  //}
  //_mtop = GetParam("mtp");
  
  // !!!!
  //for(map<string, double* >::iterator it = _xfitter_pars.begin(); it != _xfitter_pars.end(); it++)
  //{
  //  printf("_xfitter_pars[%s] = %f\n", it->first.c_str(), *it->second);
  //}

  return 0;
}

void ReactionHathor::setDatasetParameters(int dataSetID, map<std::string, std::string> pars, map<std::string, double> dsPars)
{
  // check if dataset with provided ID already exists
  if(_hathorArray.find(dataSetID) != _hathorArray.end())
  {
    char str[256];
    sprintf(str, "F: dataset with id = %d already exists", dataSetID);
    hf_errlog_(17080701, str, strlen(str));
  }

  // instantiate Hathor
  Hathor* hathor = new Hathor(*_pdf);
  //Hathor* hathor = new Hathor();

  // set collision type
  // read ppbar (0 for pp, 1 for ppbar) from provided dataset parameters
  // if not specified assume it is 0 (pp collisions)
  // here local value is preferred over global one (to allow different data sets for different collision types)
  int ppbar = false;
  if(pars.find("ppbar") != pars.end())
    ppbar = atoi(pars.find("ppbar")->second.c_str());
  else if(checkParam("ppbar"))
    ppbar = GetParamI("ppbar");
  if(ppbar)
    hathor->setColliderType(Hathor::PPBAR);
  else
    hathor->setColliderType(Hathor::PP);

  // read centre-of-mass energy from provided dataset parameters (must be provided)
  // here local value is preferred over global one (to allow different data sets with difference centre-of-mass energies)
  double sqrtS = (pars.find("SqrtS") != pars.end()) ? atof(pars.find("SqrtS")->second.c_str()) : GetParam("SqrtS");
  if(sqrtS == 0.0)
    hf_errlog(17080702, "F: no SqrtS for dataset with id = " + std::to_string(dataSetID));
  hathor->setSqrtShad(sqrtS);

  // set mass
  // here local value is preferred over global one (to allow calculations with several mass values, e.g. for translation into MSbar mass scheme)
  // the value may change further in iterations, in this case store NULL pointer (will be updated at each iteration and treated as global value)
  // TODO: fix memory leak
  _mtopPerInstance[dataSetID] = (pars.find("mtp") != pars.end()) ? new double(atof(pars.find("mtp")->second.c_str())) : NULL;

  // set renorm. scale
  _mrPerInstance[dataSetID] = _mtopPerInstance[dataSetID];
  if(_mtopPerInstance[dataSetID] && checkParam("muR"))
    *_mrPerInstance[dataSetID] *= GetParam("muR");

  // set fact. scale
  _mfPerInstance[dataSetID] = _mtopPerInstance[dataSetID];
  if(_mtopPerInstance[dataSetID] && checkParam("muF"))
    *_mfPerInstance[dataSetID] *= GetParam("muF");

  // set perturbative order
  // here local value is preferred over global one (to allow LO and NLO calculations in one run, e.g. for translation into MSbar mass scheme)
  std::string schemeName = (pars.find("Order") != pars.end()) ? pars.find("Order")->second : GetParamS("Order");
  int scheme = Hathor::LO;
  if(schemeName == "NLO")
    scheme = scheme | Hathor::NLO;
  else if(schemeName == "NNLO")
    scheme = scheme | Hathor::NNLO;
  // set mass scheme (default is pole mass scheme)
  // here local value is preferred over global one
  int msMass = 0;
  if(pars.find("MS_MASS") != pars.end())
    msMass = atoi(pars.find("MS_MASS")->second.c_str());
  else if(checkParam("MS_MASS"))
    msMass = GetParamI("MS_MASS");
  if(msMass)
    scheme = scheme | Hathor::MS_MASS;
  hathor->setScheme(scheme);

  // set precision level
  // read precision level from provided dataset parameters
  // if not specified set to default 2 -> Hathor::MEDIUM
  int precisionLevel = 2;
  if(checkParam("precisionLevel"))
    precisionLevel = GetParamI("precisionLevel");
  else if(pars.find("precisionLevel") != pars.end())
    precisionLevel = atoi(pars.find("precisionLevel")->second.c_str());
  precisionLevel = std::pow(10, 2 + precisionLevel);
  // check that this setting is allowed
  // see in AbstractHathor.h:
  //   enum ACCURACY { LOW=1000, MEDIUM=10000, HIGH=100000 };
  // and
  // precisionLevel = 1 -> Hathor::LOW
  // precisionLevel = 2 -> Hathor::MEDIUM
  // precisionLevel = 3 -> Hathor::HIGH
  if(precisionLevel !=  Hathor::LOW && precisionLevel !=  Hathor::MEDIUM && precisionLevel !=  Hathor::HIGH)
    hf_errlog(17081102, "F: provided precision level = " + std::to_string(precisionLevel) + " not supported by Hathor");
  hathor->setPrecision(precisionLevel);
  
  std::cout << " Hathor will use for this instance (" + std::to_string(dataSetID) + "):" << std::endl;
  double mt = _mtopPerInstance[dataSetID] ? *_mtopPerInstance[dataSetID] : GetParam("mtp");
  std::cout << " mtop = " << mt << "[GeV] " << std::endl;
  std::cout << " renorm. scale = " << (_mrPerInstance[dataSetID] ? *_mrPerInstance[dataSetID] : (mt * GetParam("muR"))) << "[GeV] " << std::endl;
  std::cout << " factor. scale = " << (_mfPerInstance[dataSetID] ? *_mfPerInstance[dataSetID] : (mt * GetParam("muF"))) << "[GeV] " << std::endl;
  std::cout << " SqrtS = " << sqrtS << std::endl;
  std::cout << " scheme: " << scheme << std::endl;
  std::cout << " precisionLevel: " << precisionLevel << std::endl;
  std::cout << std::endl;

  // done
  hathor->PrintOptions();
  _hathorArray[dataSetID] = hathor;
}

// Main function to compute results at an iteration
int ReactionHathor::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
  _pdf->IsValid = true;
  rlxd_reset(_rndStore);

  Hathor* hathor = _hathorArray.at(dataSetID);
  //hathor->getXsection(_mtop, _mr, _mf);
  double mt = _mtopPerInstance[dataSetID] ? *_mtopPerInstance[dataSetID] : GetParam("mtp");
  double mr = _mrPerInstance[dataSetID] ? *_mrPerInstance[dataSetID] : (mt * GetParam("muR"));
  double mf = _mfPerInstance[dataSetID] ? *_mfPerInstance[dataSetID] : (mt * GetParam("muF"));
  hathor->getXsection(mt, mr, mf);
  double dum = 0.0;
  double xsec = 0.0;
  hathor->getResult(0, xsec, dum);
  printf("mt,mr,mf,xsec: %f %f %f %f\n", mt, mr, mf, xsec);
  val = xsec;
  //printf("VAL ************ %f\n", val[0]);

  return 0;
}
