
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
  _scheme = -1;
  _mtop = -1.0;
  _mr = -1.0;
  _mf = -1.0;
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

int ReactionHathor::atStart(const string &s)
{
  // PDFs for Hathor
  _pdf = new HathorPdfxFitter(this);

  // random number generator
  rlxd_init(1, 1);
  int nRnd = rlxd_size();
  //std::cout << " Size of random number array = " << nRnd << "\n";
  _rndStore = new int[nRnd];
  rlxd_get(_rndStore);

  // scheme (perturbative order and pole/MSbar mass treatment)
  const string order = GetParamS("Order");
  const int  pertubOrder = OrderMap(order);
  _scheme = Hathor::LO;
  if(pertubOrder > 1)
    _scheme = _scheme | Hathor::NLO;
  if(pertubOrder > 2)
    _scheme = _scheme | Hathor::NNLO;
  int msMass = 0; // pole mass by default
  if(checkParam("MS_MASS"))
    msMass = GetParamI("MS_MASS");
  if(msMass)
    _scheme = _scheme | Hathor::MS_MASS;

  // top quark mass
  std::string mtopName = "mtp";// shouldn't we distinguish somehow between pole and running masses?
  if(!checkParam(mtopName))
  {
    std::string str = "F: no top quark mass (\"" + mtopName + "\" parameter) for Hathor";
    hf_errlog_(17081101, str.c_str(), strlen(str.c_str()));
  }
  _mtop = GetParam("mtp");

  // renorm. scale
  _mr = _mtop;
  if(checkParam("muR"))
    _mr *= GetParam("muR");

  // fact. scale
  _mf = _mtop;
  if(checkParam("muF"))
    _mf *= GetParam("muF");

  std::cout << " Hathor will use:";
  std::cout << " mtop = " << _mtop << "[GeV] ";
  std::cout << " renorm. scale = " << _mr << "[GeV] ";
  std::cout << " fact. scale = " << _mf << "[GeV]";
  std::cout << std::endl;

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

  // read centre-of-mass energy from provided dataset parameters
  // (must be provided)
  auto it = pars.find("SqrtS");
  if(it == pars.end())
  {
    char str[256];
    sprintf(str, "F: no SqrtS for dataset with id = %d", dataSetID);
    hf_errlog_(17080702, str, strlen(str));
  }
  double sqrtS = atof(it->second.c_str());

  // read precision level from provided dataset parameters
  // if not specified set to default 2 -> Hathor::MEDIUM
  int precisionLevel = Hathor::MEDIUM;
  it = pars.find("precisionLevel");
  if(it != pars.end())
  {
    precisionLevel = std::pow(10, 2 + atoi(it->second.c_str()));
    // check that this setting is allowed
    // see in AbstractHathor.h:
    //   enum ACCURACY { LOW=1000, MEDIUM=10000, HIGH=100000 };
    // and
    // precisionLevel = 1 -> Hathor::LOW
    // precisionLevel = 2 -> Hathor::MEDIUM
    // precisionLevel = 3 -> Hathor::HIGH
    if(precisionLevel !=  Hathor::LOW && precisionLevel !=  Hathor::MEDIUM && precisionLevel !=  Hathor::HIGH)
    {
      char str[256];
      sprintf(str, "F: provided precision level = %d not supported by Hathor", precisionLevel);
      hf_errlog_(17081102, str, strlen(str));
    }
  }

  // read ppbar from provided dataset parameters
  // if not specified assume it is false (pp collisions)
  int ppbar = false;
  it = pars.find("ppbar");
  if(it != pars.end())
  {
    ppbar = atoi(it->second.c_str());
    if(ppbar !=  0 && ppbar != 1)
    {
      char str[256];
      sprintf(str, "F: provided ppbar = %d not recognised (must be 0 or 1)", ppbar);
      hf_errlog_(17081103, str, strlen(str));
    }
  }

  // instantiate Hathor
  Hathor* hathor = new Hathor(*_pdf);
  //Hathor* hathor = new Hathor();

  // set collision type
  if(ppbar)
    hathor->setColliderType(Hathor::PPBAR);
  else
    hathor->setColliderType(Hathor::PP);

  // set centre-of-mass energy
  hathor->setSqrtShad(sqrtS);

  // set scheme
  hathor->setScheme(_scheme);

  // set precision level
  hathor->setPrecision(precisionLevel);

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
  hathor->getXsection(_mtop, _mr, _mf);
  double dum = 0.0;
  val[0] = 0.0;
  hathor->getResult(0, val[0], dum);
  //printf("VAL ************ %f\n", val[0]);

  return 0;
}

