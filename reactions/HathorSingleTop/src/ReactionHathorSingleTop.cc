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
  _scheme = -1;
  _mtop = -1.0;
  _mr = -1.0;
  _mf = -1.0;
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


// Initialize at the start of the computation
int ReactionHathorSingleTop::initAtStart(const string &s)
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
  std::string mtopName = "mtp";
  if(!checkParam(mtopName))
  {
    std::string str = "F: no top quark mass (\"" + mtopName + "\" parameter) for Hathor";
    hf_errlog_(18081701, str.c_str(), strlen(str.c_str()));
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


void ReactionHathorSingleTop::setDatasetParameters(int dataSetID, map<std::string, std::string> pars, map<std::string, double> dsPars)
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
    hf_errlog_(18081702, str, strlen(str));
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
      hf_errlog_(18081702, str, strlen(str));
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

  // read topquark from provided dataset parameters
  // if not specified assume it is false (topquark collisions)
  int antitopquark = false;
  it = pars.find("antitopquark");
  if(it != pars.end())
  {
    antitopquark = atoi(it->second.c_str());
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
  std::cout << "ReactionHathorSingleTop: TOPQUARK/ANTITOPQUARK se to " << antitopquark << std::endl;
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

  // set scheme
  hathor->setScheme(_scheme);
  std::cout << "ReactionHathorSingleTop: Setting the scheme" << std::endl;
  // set precision level
  hathor->setPrecision(precisionLevel);

  // done
  hathor->PrintOptions();
  _hathorArray[dataSetID] = hathor;
}

// Main function to compute results at an iteration
int ReactionHathorSingleTop::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
  _pdf->IsValid = true;
  rlxd_reset(_rndStore);

  HathorSgTopT* hathor = _hathorArray.at(dataSetID);

  int msMass = 0; // pole mass by default
  if(checkParam("MS_MASS"))
    msMass = GetParamI("MS_MASS");
  if(msMass){

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
 
    aspi = hathor->getAlphas(_mr)/(pi);
    dmtms = _mtop/100.;

    // decoupling coefficients
    nfl = 5.;

    Lrbar = 0.;
    d1dec = ( 4./3. + Lrbar );
    d2dec = ( 307./32. + 2.*z2 + 2./3.*z2*ln2 - z3/6.
             + 509./72.*Lrbar + 47./24.*pow(Lrbar,2)
             - nfl*(71./144. + z2/3. + 13./36.*Lrbar + pow(Lrbar,2)/12.) );
   
    
    _scheme = Hathor::LO;
    hathor->setScheme(_scheme);

    // LO
    hathor->getXsection(_mtop,_mr,_mf);
    hathor->getResult(0,valtclo,err1,chi1);
    
    // LO derivatives
    hathor->getXsection(_mtop+dmtms,_mr,_mf);
    hathor->getResult(0,valtclop,err1,chi1);
    
    hathor->getXsection(_mtop-dmtms,_mr,_mf);
    hathor->getResult(0,valtclom,err1,chi1);

    
    _scheme = Hathor::LO | Hathor::NLO;
    hathor->setScheme(_scheme) ;
  
    // NLO
    hathor->getXsection(_mtop,_mr,_mf);
    hathor->getResult(0,valtcnlo,err1,chi1);
    
    // NLO derivatives
    hathor->getXsection(_mtop+dmtms,_mr,_mf);
    hathor->getResult(0,valtcnlop,err1,chi1);
     
    hathor->getXsection(_mtop-dmtms,_mr,_mf);
    hathor->getResult(0,valtcnlom,err1,chi1);
 
    // add things up
    crst = valtcnlo
      + aspi* d1dec*_mtop/(2.*dmtms)* (valtcnlop-valtcnlom)
      + pow(aspi,2)* d2dec*_mtop/(2.*dmtms)* (valtclop-valtclom)
      + pow(aspi*d1dec*_mtop/dmtms,2)/2.* (valtclop-2.*valtclo+valtclom);

    val[0]=crst;
 
  }
  else{
  
    hathor->getXsection(_mtop, _mr, _mf);
    double dum = 0.0;
    val[0] = 0.0;
    hathor->getResult(0, val[0], dum);
    
  }
   return 0;
}

