/*
   @file ReactionHathorSingleTop.cc
   @date 2018-07-25
   @author  AddReaction.py
   Created by  AddReaction.py on 2018-07-25
*/

#include "ReactionHathorSingleTop.h"
#include "HathorPdfxFitter.h"
#include "Hathor.h"
#include "cstring"
#include "xfitter_cpp.h"

// the class factories
extern "C" ReactionHathorSingleTop* create() {
  return new ReactionHathorSingleTop();
}

extern "C"
{
int rlxd_size(void);
void rlxd_get(int state[]);
void rlxd_reset(int state[]);
void rlxd_init(int level,int seed);
}

ReactionHathorSingleTop::ReactionHathorSingleTop()
{
  _pdf = NULL;
  _rndStore = NULL;
}

ReactionHathorSingleTop::~ReactionHathorSingleTop()
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

void ReactionHathorSingleTop::initTerm(TermData *td)
{
  ReactionTheory::initTerm(td);
  int dataSetID = td->id;
  _tdDS[dataSetID] = td;

  // check if dataset with provided ID already exists
  if(_hathorArray.find(dataSetID) != _hathorArray.end())
  {
    char str[256];
    sprintf(str, "F: dataset with id = %d already exists", dataSetID);
    hf_errlog_(19060401, str, strlen(str));
  }

  // read centre-of-mass energy from provided dataset parameters
  // (must be provided)
  double sqrtS = 0.0;
  if(td->hasParam("SqrtS"))
    sqrtS = *td->getParamD("SqrtS");
  if(sqrtS == 0.0)
  {
    char str[256];
    sprintf(str, "F: no SqrtS for dataset with id = %d", dataSetID);
    hf_errlog_(18081702, str, strlen(str));
  }

  // read precision level from provided dataset parameters
  // if not specified set to default 2 -> Hathor::MEDIUM
  int precisionLevel = Hathor::MEDIUM;
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
  {
    char str[256];
    sprintf(str, "F: provided precision level = %d not supported by Hathor", precisionLevel);
    hf_errlog_(18081702, str, strlen(str));
  }

  // read ppbar from provided dataset parameters
  // if not specified assume it is false (pp collisions)
  int ppbar = false;
  if(td->hasParam("ppbar"))
  {
    ppbar = td->getParamI("ppbar");
    if(ppbar !=  0 && ppbar != 1)
    {
      char str[256];
      sprintf(str, "F: provided ppbar = %d not recognised (must be 0 or 1)", ppbar);
      hf_errlog_(17081103, str, strlen(str));
    }
  }

  // read topquark from provided dataset parameters
  // if not specified assume it is false (topquark collisions)
  int antitopquark = false;
  if(td->hasParam("antitopquark"))
  {
    antitopquark = td->getParamI("antitopquark");
    if(antitopquark !=  0 && antitopquark != 1)
    {
      char str[256];
      sprintf(str, "F: provided antitopquark = %d not recognised (must be 0 or 1)", antitopquark);
      hf_errlog_(17081103, str, strlen(str));
    }
  }

  // instantiate Hathor
  HathorSgTopT* hathor = new HathorSgTopT(*_pdf);
  hathor->sethc2(0.38937911e9);

  // set collision type
  if(ppbar)
    hathor->setColliderType(Hathor::PPBAR);
  else
    hathor->setColliderType(Hathor::PP);

  std::cout << "ReactionHathorSingleTop: PP/PPBAR parameter set to " << ppbar << std::endl;

  // set centre-of-mass energy
  hathor->setSqrtShad(sqrtS);
  std::cout << "ReactionHathorSingleTop: center of mass energy set to " << sqrtS << std::endl;
  // choose TOPQUARK/ANTITOPQUARK
  std::cout << "ReactionHathorSingleTop: TOPQUARK/ANTITOPQUARK set to " << antitopquark << std::endl;
  if(antitopquark)
  {
    std::cout << " Antitopquark is set" << std::endl;
    hathor->setParticle(SgTop::ANTITOPQUARK);
  }
  else
  {
    std::cout << " Topquark is selected" << std::endl;
    hathor->setParticle(SgTop::TOPQUARK);
  }

  // scheme (perturbative order and pole/MSbar mass treatment)
  std::string order = td->getParamS("Order");
  _scheme[dataSetID] = Hathor::LO;
  if(order == "NLO")
    _scheme[dataSetID] = _scheme[dataSetID] | Hathor::NLO;
  if(order == "NNLO")
    _scheme[dataSetID] = _scheme[dataSetID] | Hathor::NNLO;
  int msMass = 0; // pole mass by default
  if(td->hasParam("MS_MASS"))
    msMass = td->getParamI("MS_MASS");
  if(msMass)
    _scheme[dataSetID] = _scheme[dataSetID] | Hathor::MS_MASS;
  hathor->setScheme(_scheme[dataSetID]);
  std::cout << "ReactionHathorSingleTop: Setting the scheme" << std::endl;
  // set precision level
  hathor->setPrecision(precisionLevel);

  // top quark mass
  _mtop[dataSetID] = *td->getParamD("mtp");

  // renorm. scale
  _mr[dataSetID] = _mtop[dataSetID];
  if(td->hasParam("muR"))
    _mr[dataSetID] *= *td->getParamD("muR");

  // fact. scale
  _mf[dataSetID] = _mtop[dataSetID];
  if(td->hasParam("muF"))
    _mf[dataSetID] *= *td->getParamD("muF");

  std::cout << " Hathor will use:";
  std::cout << " mtop = " << _mtop[dataSetID] << "[GeV] ";
  std::cout << " renorm. scale = " << _mr[dataSetID] << "[GeV] ";
  std::cout << " fact. scale = " << _mf[dataSetID] << "[GeV]";
  std::cout << std::endl;

  // done
  hathor->PrintOptions();
  _hathorArray[dataSetID] = hathor;
}

// Initialize at the start of the computation
void ReactionHathorSingleTop::atStart()
{
  // PDFs for Hathor
  _pdf = new HathorPdfxFitter(this);

  // random number generator
  rlxd_init(1, 1);
  int nRnd = rlxd_size();
  //std::cout << " Size of random number array = " << nRnd << "\n";
  _rndStore = new int[nRnd];
  rlxd_get(_rndStore);
}
void ReactionHathorSingleTop::compute(TermData *td, valarray<double> &val, map<string, valarray<double> > &err)
{
  td->actualizeWrappers();
  _pdf->IsValid = true;
  int dataSetID = td->id;
  rlxd_reset(_rndStore);

  HathorSgTopT* hathor = _hathorArray.at(dataSetID);

  if(_scheme[dataSetID] & Hathor::MS_MASS){

    double dmtms;
    double Lrbar,nfl;
    double d1dec,d2dec,aspi;
    double const pi = 3.141592653589793;
    double const z2 = 1.644934066848226;
    double const z3 = 1.202056903159594;
    double const ln2= 0.693147180559945;

    double crst;
    double valtclo, valtclop, valtclom, valtcnlo, valtcnlop, valtcnlom;
    double err1, chi1;
    aspi = hathor->getAlphas(_mr[dataSetID])/(pi);
    dmtms = _mtop[dataSetID]/100.;

    // decoupling coefficients
    nfl = 5.;

    Lrbar = 0.;
    d1dec = ( 4./3. + Lrbar );
    d2dec = ( 307./32. + 2.*z2 + 2./3.*z2*ln2 - z3/6.
              + 509./72.*Lrbar + 47./24.*pow(Lrbar,2)
              - nfl*(71./144. + z2/3. + 13./36.*Lrbar + pow(Lrbar,2)/12.) );


    int scheme = Hathor::LO;
    hathor->setScheme(scheme);

    // LO
    hathor->getXsection(_mtop[dataSetID],_mr[dataSetID],_mf[dataSetID]);
    hathor->getResult(0,valtclo,err1,chi1);

    // LO derivatives
    hathor->getXsection(_mtop[dataSetID]+dmtms,_mr[dataSetID],_mf[dataSetID]);
    hathor->getResult(0,valtclop,err1,chi1);

    hathor->getXsection(_mtop[dataSetID]-dmtms,_mr[dataSetID],_mf[dataSetID]);
    hathor->getResult(0,valtclom,err1,chi1);


    scheme = Hathor::LO | Hathor::NLO;
    hathor->setScheme(scheme) ;

    // NLO
    hathor->getXsection(_mtop[dataSetID],_mr[dataSetID],_mf[dataSetID]);
    hathor->getResult(0,valtcnlo,err1,chi1);

    // NLO derivatives
    hathor->getXsection(_mtop[dataSetID]+dmtms,_mr[dataSetID],_mf[dataSetID]);
    hathor->getResult(0,valtcnlop,err1,chi1);

    hathor->getXsection(_mtop[dataSetID]-dmtms,_mr[dataSetID],_mf[dataSetID]);
    hathor->getResult(0,valtcnlom,err1,chi1);

    // add things up
    crst = valtcnlo
           + aspi* d1dec*_mtop[dataSetID]/(2.*dmtms)* (valtcnlop-valtcnlom)
           + pow(aspi,2)* d2dec*_mtop[dataSetID]/(2.*dmtms)* (valtclop-valtclom)
           + pow(aspi*d1dec*_mtop[dataSetID]/dmtms,2)/2.* (valtclop-2.*valtclo+valtclom);

    val[0]=crst;
    // SZ 05.05.2019 it seems that order (LO or NLO) is not treated properly above

  }
  else{

    hathor->getXsection(_mtop[dataSetID], _mr[dataSetID], _mf[dataSetID]);
    double dum = 0.0;
    val[0] = 0.0;
    hathor->getResult(0, val[0], dum);

  }
}

