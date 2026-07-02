/*
   @file ReactionHathorSingleTop.cc
   @date 2018-07-25
   @author  AddReaction.py
   Created by  AddReaction.py on 2018-07-25
*/

#include "ReactionHathorSingleTop.h"
#include "HathorPdfxFitter.h"
#include "Hathor.h"
#include "HathorGenericIntegrator.h"
#include "cstring"
#include "xfitter_cpp.h"
#include "xfitter_steer.h"
#include <unistd.h>

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
    _hathorS = NULL;
    _hathorT = NULL;
    _hathorWt = NULL;
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
    //for(auto item : _hathorTArray)
    //  if(item.second)
    //    delete item.second;
    //for(auto item : _hathorSArray)
    //  if(item.second)
    //    delete item.second;
    //for(auto item : _hathorWtArray)
    //  if(item.second)
    //    delete item.second;
}

// Compute coefficients for transforming to arbitrary alpha_s(mu) via the eq.
//   as(m_MSbar)=as(mu)(1 +  as(mu_r)*(4*pi^2)*Lrbar*bar0
//                        + as^2(mu_r)*(4*pi^2)^2*( Lrbar*bar1+Lrbar^2*bar0^2))
//Output: a vector of correction factors for LO, NLO and NNLO terms in cs.
//N.B. LO terms in cross-section have as power 2. This does not return N=1 ATM 
//but if need be, they can be computed as:
//  order > LO terms:
//    asFactor = 1.
//    asFactor += asNEW*4*pi*Lmu*bar0;
//  order > NLO extra term:
//    asFactor += pow((asNEW*4*pi),2)*( Lmu*bar1 + pow(Lmu*bar0,2) );
//  Common factor for all orders:
//    asFactor *= asNEW/asOLD;
//vector<double> ReactionHathorSingleTop::asFactors(SgTop* XS, 
vector<double> ReactionHathorSingleTop::asFactors(int msMass, double nfl, int orderI, IHathorGenericIntegrator* XS, 
                                                  double muOLD, double muNEW)
{
    vector<double> ret;   //n:th component will be the factor for as^(n+2)
    for (int n=0; n!=3; ++n) ret.push_back(1.);  //Init

    if (msMass == 0) return ret;  //No factors needed in pole scheme

    double Lmu   = log(pow(muNEW/muOLD,2));
    double asOLD = XS->getAlphas(muOLD);
    double asNEW = XS->getAlphas(muNEW);

    //Set alpha_S beta coef.s, needs orderI and nfl (#active flavors)
    double const pi = 3.141592653589793;
    double beta0  = 11. -  2.*nfl/3.;
    double beta1 = orderI > 0 ? 102. - 38.*nfl/3. : 0.;
    double bar0   = beta0/pow(4.*pi,2);
    double bar1   = beta1/pow(4.*pi,4);
    
    if (orderI != 0) {                           //O(as^3)
        ret[0] += 8.*pi*asNEW*Lmu*bar0;
        if (orderI > 1) {                        //O(as^4)
            ret[0] +=  16.*pow(pi*asNEW*Lmu*bar0,2)
                     + 32.*pow(pi*asNEW,2)*( Lmu*bar1 + pow(Lmu*bar0,2) );
            ret[1]  = 1. + 12.*pi*asNEW*Lmu*bar0;
            //All additions to ret[2] would be O(as^5)
        }
    }
    for (unsigned int n=0; n!=ret.size(); ++n) ret[n] *= pow(asNEW/asOLD, n+2.);
    
    return ret;
}

void ReactionHathorSingleTop::initTerm(TermData *td)
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
  _alphaem = td->getParamD("alphaem");
  _sin2thW = td->getParamD("sin2thW");
  _ckm[0][0] = *td->getParamD("Vud");
  _ckm[0][1] = *td->getParamD("Vus");
  _ckm[0][2] = *td->getParamD("Vub");
  _ckm[1][0] = *td->getParamD("Vcd");
  _ckm[1][1] = *td->getParamD("Vcs");
  _ckm[1][2] = *td->getParamD("Vcb");
  _ckm[2][0] = *td->getParamD("Vtd");
  _ckm[2][1] = *td->getParamD("Vts");
  _ckm[2][2] = *td->getParamD("Vtb");
  _mtopPerInstance[dataSetID] = td->getParamD("mtp");
  _mfPerInstance[dataSetID] = *td->getParamD("muR");
  _mrPerInstance[dataSetID] = *td->getParamD("muF");
  _orderPerInstance[dataSetID] = td->getParamS("Order");
  _MsbarPerInstance[dataSetID] = td->getParamI("MS_MASS");
  _precisionLevelPerInstance[dataSetID] = td->getParamI("precisionLevel");
  _sqrtSPerInstance[dataSetID] = *td->getParamD("SqrtS");
  _ppbarPerInstance[dataSetID] = td->getParamI("ppbar");

    // read enabled processes. By default, enable everything
    _tchannel[dataSetID]=1;
    _schannel[dataSetID]=1;
    _Wtchannel[dataSetID]=1;
    if(td->hasParam("tSgTop") ) _tchannel[dataSetID]  = td->getParamI("tSgTop");
    if(td->hasParam("sSgTop") ) _schannel[dataSetID]  = td->getParamI("sSgTop");
    if(td->hasParam("WtSgTop")) _Wtchannel[dataSetID] = td->getParamI("WtSgTop");
    if (_tchannel[dataSetID]!=1 && _schannel[dataSetID]!=1 && _Wtchannel[dataSetID]!=1) {
        hf_errlog(21121001,"F: ERROR all channels disabled in ReactionHathorSingleTop");
    } else {
        if (_tchannel[dataSetID] ==0) hf_errlog(21121002,"I: Disabled t-channel processes in ReactionHathorSingleTop");
        if (_schannel[dataSetID] ==0) hf_errlog(21121003,"I: Disabled s-channel processes in ReactionHathorSingleTop");
        if (_Wtchannel[dataSetID]==0) hf_errlog(21121004,"I: Disabled W+t final state processes in ReactionHathorSingleTop");
    }
    _antiquark[dataSetID] = 0;
    if(td->hasParam("antitopquark")) {
        _antiquark[dataSetID] = td->getParamI("antitopquark");
        if(_antiquark[dataSetID] !=  0 && _antiquark[dataSetID] != 1) {
            hf_errlog(17081103, "F: provided antitopquark = " + std::to_string(_antiquark[dataSetID]) + " not recognised (must be 0 or 1)");
        }
    }

  // avoid calculating the same cross section many times
  std::string calc_name = "mtp" + std::to_string((unsigned long long)(void**)_mtopPerInstance[dataSetID]);
  calc_name += std::string("_") + "muF" + std::to_string(_mfPerInstance[dataSetID]);
  calc_name += std::string("_") + "muR" + std::to_string(_mrPerInstance[dataSetID]);
  calc_name += std::string("_") + "Order" + _orderPerInstance[dataSetID];
  calc_name += std::string("_") + "MS_MASS" + std::to_string(_MsbarPerInstance[dataSetID]);
  calc_name += std::string("_") + "precisionLevel" + std::to_string(_precisionLevelPerInstance[dataSetID]);
  calc_name += std::string("_") + "SqrtS" + std::to_string(_sqrtSPerInstance[dataSetID]);
  calc_name += std::string("_") + "ppbar" + std::to_string(_ppbarPerInstance[dataSetID]);
  calc_name += std::string("_") + "tchannel" + std::to_string(_tchannel[dataSetID]);
  calc_name += std::string("_") + "schannel" + std::to_string(_schannel[dataSetID]);
  calc_name += std::string("_") + "Wtchannel" + std::to_string(_Wtchannel[dataSetID]);
  calc_name += std::string("_") + "antitopquark" + std::to_string(_antiquark[dataSetID]);
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
  //printf("adding calc_name = %s td = %p\n", calc_name.c_str(), td);
  _convolved[calc_name].second.push_back(td);
}

// Initialize at the start of the computation
void ReactionHathorSingleTop::atStart()
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

    _hathorS = new HathorGenericIntegrator<HathorSgTopS>(*_pdf);
    _hathorT = new HathorGenericIntegrator<HathorSgTopT>(*_pdf);
    _hathorWt = new HathorGenericIntegrator<HathorSgTopWt>(*_pdf);
}
void ReactionHathorSingleTop::compute(TermData *td, valarray<double> &val, map<string, valarray<double> > &err)
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

void ReactionHathorSingleTop::atIteration() {
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

    // perturbative order and mass scheme
    int scheme = Hathor::LO;
    int orderI = 0;
    double kfactors_nnlo_tch = 1.0;
    if (_orderPerInstance[dataSetID] == "NLO") {
        scheme = scheme | Hathor::NLO;  
        orderI = 1;                                             
    } else if (_orderPerInstance[dataSetID] == "NNLO") {
        hf_errlog(21121005,"W: Standard Hathor-2.0 has no NNLO single top processes. ReactionHathorSingleTop reverts to NLO.");
        scheme = scheme | Hathor::NLO;  
        orderI = 1;                                             
        // NNLO/NLO K-factor as given in ABMP16 paper arXiv:1701.05838 with refs. to arXiv:1404.7116 and arXiv:1608.05212
        kfactors_nnlo_tch = 0.984;
        /* NNLO not implemented in Hathor-2.0. If updated, remove the above 3 
            *lines and uncomment the below */
        //scheme = scheme | Hathor::NLO  | Hathor::NNLO;
        //orderI = 2;                                             
    } else if (_orderPerInstance[dataSetID] != "LO") {
        //std::cout << " ReactionHathorSingleTop: perturbative order "
        //                        << order 
        //                        <<  " not supported. Defaulting to NLO."
        //                        << std::endl;
        scheme = scheme | Hathor::NLO;
        orderI = 1;                                             
    }
    scheme = scheme | Hathor::MS_MASS;

    // Hathor objects for different processes
    vector<IHathorGenericIntegrator*> hathorChannels;
    if (_tchannel[dataSetID]) {
        hathorChannels.push_back(_hathorT);
    }
    if (_schannel[dataSetID]) {
        hathorChannels.push_back(_hathorS);
    }
    if (_Wtchannel[dataSetID]) {
        hathorChannels.push_back(_hathorWt);
    }
    for (auto hathor : hathorChannels) {
        // collision type
        if(_ppbarPerInstance[dataSetID])
        hathor->setColliderType(Hathor::PPBAR);
        else
        hathor->setColliderType(Hathor::PP);

        // read centre-of-mass energy
        hathor->setSqrtShad(_sqrtSPerInstance[dataSetID]);

        // set conversion factor
        static constexpr double convFac_in = 0.389379323e9;  //HATHOR default value
        hathor->sethc2(convFac_in);

        // set EW constants
        hathor->setSwq(*_sin2thW);
        hathor->setAlpha(*_alphaem);
        hathor->setCkmMatrix(_ckm);

        // perturbative order and mass scheme
        hathor->setScheme(scheme);

        // set particle
        if(_antiquark[dataSetID]) {
            hathor->setParticle(SgTop::ANTITOPQUARK);
        }
        else {
            hathor->setParticle(SgTop::TOPQUARK);
        }

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
        hathor->setPrecision(precisionLevel);

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
            hathor->PrintOptions();
        }
    }

    xsec = 0.0;
    for (auto hathor : hathorChannels) {

        double crst=0.;  //Total inclusive cross section for one process

        if (_MsbarPerInstance[dataSetID] != 0) {
    
            double valtclo, valtclop, valtclom, valtcnlo, valtcnlop, valtcnlom, valtcnnlo;
            double err1, chi1;
            double aspi  = hathor->getAlphas(mt)/(pi);
            double dmtms = mt/100.;  //For numerical derivative w.r.t mt
    
            // decoupling coefficients
            static constexpr double nfl = 5.;
            static constexpr double Lrbar = 0.;  //The logarithm is zero when mt evaluated at mu_r = mu_m = mt (*)
            static constexpr double d1dec = ( 4./3. + Lrbar );
            static constexpr double d2dec = ( 307./32. + 2.*z2 + 2./3.*z2*ln2 - z3/6.
                           + 509./72.*Lrbar + 47./24.*pow(Lrbar,2)
                           - nfl*(71./144. + z2/3. + 13./36.*Lrbar + pow(Lrbar,2)/12.) );
    
            // Use numerical stencil for MSbar transformation
            //(*) requires mu_r to be se to mt in all getXsection calls
        
            // LO
            //hathor->setScheme(Hathor::LO);
            hathor->setScheme(scheme);
            // all order result
            hathor->getXsection(mt,mt,mf, td->getParamS("integrator"));
            hathor->getResult(0,valtclo,err1,chi1);
    
            if (orderI > 0) {
                // LO derivatives
                hathor->setScheme(Hathor::LO);
                hathor->getXsection(mt+dmtms,mt,mf, td->getParamS("integrator"));
                hathor->getResult(0,valtclop,err1,chi1);
                hathor->getXsection(mt-dmtms,mt,mf, td->getParamS("integrator"));
                hathor->getResult(0,valtclom,err1,chi1);
    
                // NLO
                //hathor->setScheme(Hathor::NLO);
                //hathor->getXsection(mt,mt,mf, td->getParamS("integrator"));
                //hathor->getResult(0,valtcnlo,err1,chi1);
            }
    
            if (orderI > 1) {
                // NLO derivatives
                hathor->setScheme(Hathor::NLO);
                hathor->getXsection(mt+dmtms,mt,mf, td->getParamS("integrator"));
                hathor->getResult(0,valtcnlop,err1,chi1);
                hathor->getXsection(mt-dmtms,mt,mf, td->getParamS("integrator"));
                hathor->getResult(0,valtcnlom,err1,chi1);
    
                // NNLO
                //hathor->setScheme(Hathor::NNLO);
                //hathor->getXsection(mt,mt,mf, td->getParamS("integrator"));
                //hathor->getResult(0,valtcnnlo,err1,chi1);
            }
    
            //Coefficients for generalizing cross-section to arbitrary alpha_s(mu_r)
            vector<double> asFac = asFactors(_MsbarPerInstance[dataSetID], nfl, orderI, hathor,mt,mr);
            if (asFac.size()!=3) {
                hf_errlog(21120901,"F: ERROR in calculating as conversion factors in ReactionHathorSingleTop.cc");
                return;
            }
            double asLO   = asFac[0];
            double asNLO  = asFac[1];
            double asNNLO = asFac[2];
        
            //Combine terms to get cross-section
            double NLOder=0., NNLOder=0.;
            crst =                   asLO*valtclo           //Common LO
                   + (orderI > 0 ?  asNLO*valtcnlo  : 0.)   //Common NLO
                   + (orderI > 1 ? asNNLO*valtcnnlo : 0.);  //Common NNLO
            //Numerical derivative contributions
            if (orderI > 0) {
                NLOder = aspi*d1dec*(mt)/(2.*dmtms)*(valtclop-valtclom);
                crst  += asNLO*NLOder;
            }
            if (orderI > 1) {
                //N.B. csNLO terms include one factor of aspi on 2nd line
                NNLOder = pow(aspi,2)*d2dec*(mt)/(2.*dmtms)*(valtclop - valtclom)                  
                         +       aspi*d1dec*(mt)/(2.*dmtms)*(valtcnlop - valtcnlom)
                         +pow(aspi*d1dec*(mt)/dmtms,2)/2.*(valtclop - 2.*valtclo + valtclom);
                crst  += asNNLO*NNLOder;            
            }
            //printf("xsec,err: %f %f [%.3f%%]\n", crst, err1, err1/crst*100.);
                       
        } else {  //POLE scheme calculated in Hathor, no ext. numerical derivatives
    
            hathor->getXsection(mt, mr, mf, td->getParamS("integrator"));
            double dum = 0.0;
            hathor->getResult(0, crst, dum);
            //printf("xsec,err: %f %f [%.3f%%]\n", crst, dum, dum/crst*100.);
    
        }

        xsec += crst;

    }
    if (_tchannel[dataSetID] == 1 && _schannel[dataSetID] == 0 && _Wtchannel[dataSetID] == 0) {
        xsec *= kfactors_nnlo_tch;
    }
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
