
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
#include "xfitter_cpp.h"
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
  // PDFs for Hathor
  _pdf = new HathorPdfxFitter(this);

    // random number generator
  rlxd_init(1, 1);
  int nRnd = rlxd_size();
  //std::cout << " Size of random number array = " << nRnd << "\n";
  _rndStore = new int[nRnd];
  rlxd_get(_rndStore);

  // instantiate one Hathor instance for all terms
  _hathor = new Hathor(*_pdf);
}

void ReactionHathor::initTerm(TermData *td)
{
  //These two lines were added by Sasha Zenaiev but they are not needed
  //ReactionTheory::initTerm(td);
  //_tdDS[dataSetID] = td;

  int dataSetID = td->id;

  // check if dataset with provided ID already exists
  if(_mtopPerInstance.find(dataSetID) != _mtopPerInstance.end() && !chi2scan_.scan)
  {
    char str[256];
    sprintf(str, "F: dataset with id = %d already exists", dataSetID);
    hf_errlog_(17080701, str, strlen(str));
  }

  // instantiate Hathor for each term
  //Hathor* hathor = new Hathor(*_pdf);
  //_hathorArray[dataSetID] = _hathor;
}

// Main function to compute results at an iteration
void ReactionHathor::compute(TermData *td, valarray<double> &val, map<string, valarray<double> > &err)
{
  td->actualizeWrappers();
  int dataSetID = td->id;
  _pdf->IsValid = true;
  rlxd_reset(_rndStore);

  //Suppress Hathor output
  // SZ 2023.11.11 freopen breaks pipe redicrection, see https://c-faq.com/stdio/undofreopen.html
  // Solution from https://stackoverflow.com/questions/1908687/how-to-redirect-the-output-back-to-the-screen-after-freopenout-txt-a-stdo
  int o;
  if (!steering_.ldebug) {
    o = dup(fileno(stdout));
    freopen("/dev/null", "a", stdout);
  }
  
  //_hathor->getXsection(_mtop, _mr, _mf);

  // set collision type
  // read ppbar (0 for pp, 1 for ppbar) from provided dataset parameters
  // if not specified assume it is 0 (pp collisions)
  // here local value is preferred over global one (to allow different data sets for different collision types)
  int ppbar = false;
  if(td->hasParam("ppbar"))
    ppbar = td->getParamI("ppbar");
  else if(td->hasParam("ppbar"))
    ppbar = td->getParamI("ppbar");
  if(ppbar)
    _hathor->setColliderType(Hathor::PPBAR);
  else
    _hathor->setColliderType(Hathor::PP);

  // read centre-of-mass energy from provided dataset parameters (must be provided)
  // here local value is preferred over global one (to allow different data sets with difference centre-of-mass energies)
  double sqrtS = 0.0;
  if(td->hasParam("SqrtS"))
    sqrtS = *td->getParamD("SqrtS");
  else
    hf_errlog(17080702, "F: no SqrtS for dataset with id = " + std::to_string(dataSetID));
  _hathor->setSqrtShad(sqrtS);

  // set conversion factor
  double convFac_in = 0.389379323e9;  //HATHOR default value
  if(td->hasParam("convFac")) convFac_in = *td->getParamD("convFac");
  _hathor->sethc2(convFac_in);
  std::cout << " ReactionHathor: hc2 set to " << convFac_in << std::endl;

  // set mass
  double mtop = 0;
  if(td->hasParam("mtp"))
    mtop = *td->getParamD("mtp");
  
  // here local value is preferred over global one (to allow calculations with several mass values, e.g. for translation into MSbar mass scheme)
  _mtopPerInstance[dataSetID] = std::shared_ptr<double>(new double(mtop));
  if(td->hasParam("mtp"))
    *_mtopPerInstance[dataSetID] = mtop;

  // set renorm. scale
  _mrPerInstance[dataSetID] = std::shared_ptr<double>(new double(*_mtopPerInstance[dataSetID]));
  if(td->hasParam("muR"))
    *_mrPerInstance[dataSetID] *= *td->getParamD("muR");

  // set fact. scale
  _mfPerInstance[dataSetID] = std::shared_ptr<double>(new double(*_mtopPerInstance[dataSetID]));
  if(td->hasParam("muF"))
    *_mfPerInstance[dataSetID] *= *td->getParamD("muF");

  // set perturbative order
  // here local value is preferred over global one (to allow LO and NLO calculations in one run, e.g. for translation into MSbar mass scheme)
  std::string schemeName = td->getParamS("Order");
  int scheme = Hathor::LO;
  if(schemeName == "NLO")
    scheme = scheme | Hathor::NLO;
  else if(schemeName == "NNLO")
    scheme = scheme | Hathor::NLO | Hathor::NNLO;
  // set mass scheme (default is pole mass scheme)
  // here local value is preferred over global one
  int msMass = 0; // pole mass by default
  if(td->hasParam("MS_MASS"))
    msMass = td->getParamI("MS_MASS");
  if(msMass)
    scheme = scheme | Hathor::MS_MASS;
  _hathor->setScheme(scheme);

  // set precision level
  // read precision level from provided dataset parameters
  // if not specified set to default 2 -> Hathor::MEDIUM
  int precisionLevel = 2;
  if(td->hasParam("precisionLevel"))
    precisionLevel = td->getParamI("precisionLevel");
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
  _hathor->setPrecision(precisionLevel);

  //Resume standard output
  if (!steering_.ldebug)
  {
    //freopen ("/dev/tty", "a", stdout);
    dup2(o,fileno(stdout));
    close(o);
  }

  if (steering_.ldebug)
    {
      std::cout << " Hathor will use for this instance (" + std::to_string(dataSetID) + "):" << std::endl;
      double mt = *_mtopPerInstance[dataSetID];

      std::cout << " mtop = " << mt << "[GeV] " << std::endl;
      std::cout << " renorm. scale = " << *_mrPerInstance[dataSetID] << "[GeV] " << std::endl;
      std::cout << " factor. scale = " << *_mfPerInstance[dataSetID] << "[GeV] " << std::endl;
      std::cout << " SqrtS = " << sqrtS << std::endl;
      std::cout << " scheme: " << scheme << std::endl;
      std::cout << " precisionLevel: " << precisionLevel << std::endl;
      std::cout << std::endl;

      // done
      _hathor->PrintOptions();
    }

  double mt = _mtopPerInstance[dataSetID] ? (*_mtopPerInstance[dataSetID]) : mtop;
  double mr = *_mrPerInstance[dataSetID];
  double mf = *_mfPerInstance[dataSetID];
  _hathor->getXsection(mt, mr, mf);
  double dum = 0.0;
  double xsec = 0.0;
  _hathor->getResult(0, xsec, dum);
  //printf("mt,mr,mf,xsec: %f %f %f %f\n", mt, mr, mf, xsec);
  val = xsec;
  //printf("VAL ************ %f\n", val[0]);
}
