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
    _tdDS[dataSetID] = td;

    // check if dataset with provided ID already exists
    //if(_hathorTArray.find( dataSetID) != _hathorTArray.end() ||
    //   _hathorSArray.find( dataSetID) != _hathorSArray.end() ||
    //   _hathorWtArray.find(dataSetID) != _hathorWtArray.end()  )
    //{
    //    char str[256];
    //    sprintf(str, "F: dataset with id = %d already exists", dataSetID);
    //    hf_errlog(19060401, str);
    //}
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

    _hathorS = new HathorGenericIntegrator<HathorSgTopS>(*_pdf);
    _hathorT = new HathorGenericIntegrator<HathorSgTopT>(*_pdf);
    _hathorWt = new HathorGenericIntegrator<HathorSgTopWt>(*_pdf);
}
void ReactionHathorSingleTop::compute(TermData *td, valarray<double> &val, map<string, valarray<double> > &err)
{
    // avoid calculating the same cross section again
    std::string calc_name = "";
    std::vector<std::string> v_double_pars = {"SqrtS", "convFac", "mtp", "muR", "muF"};
    std::vector<std::string> v_int_pars = {"precisionLevel", "MS_MASS", "NFlavour", "ppbar", "tSgTop", "sSgTop", "WtSgTop", "antitopquark"};
    std::vector<std::string> v_str_pars = {"Order", "evolution", "evolution1", "evolution2", "integrator"};
    for (auto& p : v_double_pars) {
        if(td->hasParam(p)) {
            calc_name += "_" + p + std::to_string(*td->getParamD(p));
        }
    }
    for (auto& p : v_int_pars) {
        if(td->hasParam(p)) {
            calc_name += "_" + p + std::to_string(td->getParamI(p));
        }
    }
    for (auto& p : v_str_pars) {
        if(td->hasParam(p)) {
            calc_name += "_" + p + td->getParamS(p);
        }
    }
    if (_convolved.find(calc_name) != _convolved.end()) {
        val = _convolved[calc_name];
        return;
    }

    td->actualizeWrappers();
    _pdf->IsValid = true;
    int dataSetID = td->id;
    rlxd_reset(_rndStore);

    //Suppress Hathor output
    int o;
    if (!steering_.ldebug) {
        o = dup(fileno(stdout));
        if (freopen("/dev/null", "a", stdout)) {
        }
    }
    
    // read centre-of-mass energy from provided dataset parameters
    // (must be provided)
    double sqrtS = 0.0;
    if(td->hasParam("SqrtS")) sqrtS = *td->getParamD("SqrtS");
    if(sqrtS == 0.0) {
        char str[256];
        sprintf(str, "F: no SqrtS for dataset with id = %d", dataSetID);
        hf_errlog(18081702, str);
    }

    // read precision level from provided dataset parameters
    // if not specified set to default 2 -> Hathor::MEDIUM
    int precisionLevel = Hathor::MEDIUM;
    if(td->hasParam("precisionLevel")) {
        precisionLevel = td->getParamI("precisionLevel");
        precisionLevel = std::pow(10, 2 + precisionLevel);
    }
    // check that this setting is allowed
    // see in AbstractHathor.h:
    //   enum ACCURACY { LOW=1000, MEDIUM=10000, HIGH=100000 };
    // and
    // precisionLevel = 1 -> Hathor::LOW
    // precisionLevel = 2 -> Hathor::MEDIUM
    // precisionLevel = 3 -> Hathor::HIGH
    if (precisionLevel !=  Hathor::LOW    && 
        precisionLevel !=  Hathor::MEDIUM && 
        precisionLevel !=  Hathor::HIGH     )
    {
        char str[256];
        sprintf(str, "F: provided precision level = %d not supported by Hathor", precisionLevel);
        hf_errlog(18081702, str);
    }

    // read ppbar from provided dataset parameters
    // if not specified assume it is false (pp collisions)
    int ppbar = false;
    if(td->hasParam("ppbar")) {
        ppbar = td->getParamI("ppbar");
        if(ppbar !=  0 && ppbar != 1) {
            char str[256];
            sprintf(str, "F: provided ppbar = %d not recognised (must be 0 or 1)", ppbar);
            hf_errlog(17081103, str);
        }
    }

    // read topquark from provided dataset parameters
    // if not specified assume it is false (topquark collisions)
    int antitopquark = false;
    if(td->hasParam("antitopquark")) {
        antitopquark = td->getParamI("antitopquark");
        if(antitopquark !=  0 && antitopquark != 1) {
            char str[256];
            sprintf(str, "F: provided antitopquark = %d not recognised (must be 0 or 1)", antitopquark);
            hf_errlog(17081103, str);
        }
    }

    // read enabled processes. By default, enable everything
    int tchannel=1;
    int schannel=1;
    int Wtchannel=1;
    if(td->hasParam("tSgTop") ) tchannel  = td->getParamI("tSgTop");
    if(td->hasParam("sSgTop") ) schannel  = td->getParamI("sSgTop");
    if(td->hasParam("WtSgTop")) Wtchannel = td->getParamI("WtSgTop");
    if (tchannel!=1 && schannel!=1 && Wtchannel!=1) {
        hf_errlog(21121001,"F: ERROR all channels disabled in ReactionHathorSingleTop");
	} else {
        if (tchannel ==0) hf_errlog(21121002,"I: Disabled t-channel processes in ReactionHathorSingleTop");
        if (schannel ==0) hf_errlog(21121003,"I: Disabled s-channel processes in ReactionHathorSingleTop");
        if (Wtchannel==0) hf_errlog(21121004,"I: Disabled W+t final state processes in ReactionHathorSingleTop");
	}
    
    // Hathor objects for different processes
    vector<IHathorGenericIntegrator*> hathorChannels;
    if (tchannel) {
        hathorChannels.push_back(_hathorT);
    }
    if (schannel) {
        hathorChannels.push_back(_hathorS);
    }
    if (Wtchannel) {
        hathorChannels.push_back(_hathorWt);
    }
    
    int msMass = 0; // pole mass by default

    // top quark mass
    const double* mtop = td->getParamD("mtp");

    // renorm. scale
    const double* mr = td->getParamD("muR");

    // fact. scale
    const double* mf = td->getParamD("muF");

    std::string order = td->getParamS("Order");
    int scheme = Hathor::LO;  //POLE uses this
    int orderI = 0;  //MSBAR uses this
    double kfactors_nnlo_tch = 1.0;
    
    for (auto hathor : hathorChannels) {
    
        // set collision type
        if (ppbar) hathor->setColliderType(Hathor::PPBAR);
        else       hathor->setColliderType(Hathor::PP);
    
        if (init1) std::cout << " ReactionHathorSingleTop: PP/PPBAR parameter set to "
                             << ppbar << std::endl;

        // set conversion factor
        double convFac_in = 0.38937911e9;  //MCFM value by default
        if(td->hasParam("convFac")) convFac_in = *td->getParamD("convFac");
        hathor->sethc2(convFac_in);
        std::cout << " ReactionHathorSingleTop: hc2 set to "
                  << convFac_in << std::endl;

        // set EW parameters
        double sin2thW_in = 0.2228972;  //HATHOR default value
        if (td->hasParam("sin2thW")) {
            sin2thW_in = *td->getParamD("sin2thW");
            hathor->setSwq(sin2thW_in);
            std::cout << " ReactionHathorSingleTop: Swq set to "
                      << sin2thW_in << std::endl;
        }
        double alphaem_in = 1. / 132.2332298;  //HATHOR default value
        if (td->hasParam("alphaem")) {
            alphaem_in = *td->getParamD("alphaem");
            hathor->setAlpha(alphaem_in);
            std::cout << " ReactionHathorSingleTop: Alpha set to "
                      << alphaem_in << std::endl;
        }
        
        // set CKM matrix
        double ckm[3][3];
        hathor->getCkmMatrix(ckm);
        if (td->hasParam("Vud")) ckm[0][0] = *td->getParamD("Vud");
        if (td->hasParam("Vus")) ckm[0][1] = *td->getParamD("Vus");
        if (td->hasParam("Vub")) ckm[0][2] = *td->getParamD("Vub");
        if (td->hasParam("Vcd")) ckm[1][0] = *td->getParamD("Vcd");
        if (td->hasParam("Vcs")) ckm[1][1] = *td->getParamD("Vcs");
        if (td->hasParam("Vcb")) ckm[1][2] = *td->getParamD("Vcb");
        if (td->hasParam("Vtd")) ckm[2][0] = *td->getParamD("Vtd");
        if (td->hasParam("Vts")) ckm[2][1] = *td->getParamD("Vts");
        if (td->hasParam("Vtb")) ckm[2][2] = *td->getParamD("Vtb");
        hathor->setCkmMatrix(ckm);
        hathor->PrintCkmMatrix();
        
        // set centre-of-mass energy
        hathor->setSqrtShad(sqrtS);
        if (init1) std::cout << " ReactionHathorSingleTop: center of mass energy set to "
                             << sqrtS << std::endl;
        // choose TOPQUARK/ANTITOPQUARK
        if (init1) std::cout << " ReactionHathorSingleTop: TOPQUARK/ANTITOPQUARK set to "
                             << antitopquark << std::endl;
        if(antitopquark) {
            if (init1) std::cout << " Antitopquark is set" << std::endl;
            hathor->setParticle(SgTop::ANTITOPQUARK);
        } else {
            if (init1) std::cout << " Topquark is selected" << std::endl;
            hathor->setParticle(SgTop::TOPQUARK);
        }
    
        // scheme (perturbative order and pole/MSbar mass treatment)
        if (order == "NLO") {
            scheme = scheme | Hathor::NLO;  
            orderI = 1;                                             
        } else if (order == "NNLO") {
            hf_errlog(21121005,"W: Standard Hathor-2.0 has no NNLO single top processes. ReactionHathorSingleTop reverts to NLO.");
            scheme = scheme | Hathor::NLO;  
            orderI = 1;                                             
            // NNLO/NLO K-factor as given in ABMP16 paper arXiv:1701.05838 with refs. to arXiv:1404.7116 and arXiv:1608.05212
            kfactors_nnlo_tch = 0.984;
			/* NNLO not implemented in Hathor-2.0. If updated, remove the above 3 
			 *lines and uncomment the below */
            //scheme = scheme | Hathor::NLO  | Hathor::NNLO;
            //orderI = 2;                                             
        } else if (order != "LO") {
            if (init1) std::cout << " ReactionHathorSingleTop: perturbative order "
                                 << order 
                                 <<  " not supported. Defaulting to NLO."
                                 << std::endl;
            scheme = scheme | Hathor::NLO;
            orderI = 1;                                             
        }
        if(td->hasParam("MS_MASS")) msMass = td->getParamI("MS_MASS");
        if(msMass) scheme = scheme | Hathor::MS_MASS;
        hathor->setScheme(scheme);
        if (init1) std::cout << "ReactionHathorSingleTop: Setting the scheme"
                             << std::endl;
    
        // set precision level
        hathor->setPrecision(precisionLevel);
        
        if (init1) {
            std::cout << " Hathor will use:";
            std::cout << " mtop = " << *mtop << "[GeV] ";
            std::cout << " renorm. scale = " << *mr * (*mtop) << "[GeV] ";
            std::cout << " fact. scale = " << *mf * (*mtop) << "[GeV]";
            std::cout << std::endl;
        } 
    
    }

    if (init1) {
        if (tchannel) {
            _hathorT->PrintOptions(); 
        }
        if (schannel) {
            _hathorS->PrintOptions(); 
        }
        if (Wtchannel) {
            _hathorWt->PrintOptions();
        }
    }

    //Resume standard output
    if (!steering_.ldebug)
    {
        //freopen ("/dev/tty", "a", stdout);
        dup2(o,fileno(stdout));
        close(o);
    }

    val[0] = 0.;  //Final result will be stored here
    
    for (auto hathor : hathorChannels) {

        double crst=0.;  //Total inclusive cross section for one process

        if (msMass != 0) {
    
            double valtclo, valtclop, valtclom, valtcnlo, valtcnlop, valtcnlom, valtcnnlo;
            double err1, chi1;
            double aspi  = hathor->getAlphas(*mtop)/(pi);
            double dmtms = *mtop/100.;  //For numerical derivative w.r.t mt
    
            // decoupling coefficients
            double nfl = 5.;
            double Lrbar = 0.;  //The logarithm is zero when mt evaluated at mu_r = mu_m = mt (*)
            double d1dec = ( 4./3. + Lrbar );
            double d2dec = ( 307./32. + 2.*z2 + 2./3.*z2*ln2 - z3/6.
                           + 509./72.*Lrbar + 47./24.*pow(Lrbar,2)
                           - nfl*(71./144. + z2/3. + 13./36.*Lrbar + pow(Lrbar,2)/12.) );
    
            // Use numerical stencil for MSbar transformation
            //(*) requires mu_r to be se to mt in all getXsection calls
        
            // LO
            //hathor->setScheme(Hathor::LO);
            hathor->setScheme(scheme);
            // all order result
            hathor->getXsection(*mtop,*mtop,*mf*(*mtop), td->getParamS("integrator"));
            hathor->getResult(0,valtclo,err1,chi1);
    
            if (orderI > 0) {
                // LO derivatives
                hathor->setScheme(Hathor::LO);
                hathor->getXsection(*mtop+dmtms,*mtop,*mf*(*mtop), td->getParamS("integrator"));
                hathor->getResult(0,valtclop,err1,chi1);
                hathor->getXsection(*mtop-dmtms,*mtop,*mf*(*mtop), td->getParamS("integrator"));
                hathor->getResult(0,valtclom,err1,chi1);
    
                // NLO
                //hathor->setScheme(Hathor::NLO);
                //hathor->getXsection(*mtop,*mtop,*mf*(*mtop), td->getParamS("integrator"));
                //hathor->getResult(0,valtcnlo,err1,chi1);
            }
    
            if (orderI > 1) {
                // NLO derivatives
                hathor->setScheme(Hathor::NLO);
                hathor->getXsection(*mtop+dmtms,*mtop,*mf*(*mtop), td->getParamS("integrator"));
                hathor->getResult(0,valtcnlop,err1,chi1);
                hathor->getXsection(*mtop-dmtms,*mtop,*mf*(*mtop), td->getParamS("integrator"));
                hathor->getResult(0,valtcnlom,err1,chi1);
    
                // NNLO
                //hathor->setScheme(Hathor::NNLO);
                //hathor->getXsection(*mtop,*mtop,*mf*(*mtop), td->getParamS("integrator"));
                //hathor->getResult(0,valtcnnlo,err1,chi1);
            }
    
            //Coefficients for generalizing cross-section to arbitrary alpha_s(mu_r)
            vector<double> asFac = asFactors(msMass, nfl, orderI, hathor,*mtop,*mr*(*mtop));
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
                NLOder = aspi*d1dec*(*mtop)/(2.*dmtms)*(valtclop-valtclom);
                crst  += asNLO*NLOder;
            }
            if (orderI > 1) {
                //N.B. csNLO terms include one factor of aspi on 2nd line
                NNLOder = pow(aspi,2)*d2dec*(*mtop)/(2.*dmtms)*(valtclop - valtclom)                  
                         +       aspi*d1dec*(*mtop)/(2.*dmtms)*(valtcnlop - valtcnlom)
                         +pow(aspi*d1dec*(*mtop)/dmtms,2)/2.*(valtclop - 2.*valtclo + valtclom);
                crst  += asNNLO*NNLOder;            
            }
            //printf("xsec,err: %f %f [%.3f%%]\n", crst, err1, err1/crst*100.);
                       
        } else {  //POLE scheme calculated in Hathor, no ext. numerical derivatives
    
            hathor->getXsection(*mtop, *mr*(*mtop), *mf*(*mtop), td->getParamS("integrator"));
            double dum = 0.0;
            hathor->getResult(0, crst, dum);
            //printf("xsec,err: %f %f [%.3f%%]\n", crst, dum, dum/crst*100.);
    
        }

        val[0] += crst;

    }
    if (tchannel == 1 && schannel == 0 && Wtchannel == 0) {
        val[0] *= kfactors_nnlo_tch;
    }

    _convolved.insert(std::make_pair(calc_name, val[0]));

    init1 = false;
}

void ReactionHathorSingleTop::atIteration() {
    //printf("clearing %ld Hathor saved xsecs\n", _convolved.size());
    _convolved.clear();
}

bool ReactionHathorSingleTop::init1 = true;
