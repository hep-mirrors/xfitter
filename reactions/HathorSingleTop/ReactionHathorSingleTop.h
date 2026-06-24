#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionHathorSingleTop

  @brief A wrapper class for HathorSingleTop reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2019-06-04
// Authors: Laia Parets Peris <laia.parets.peris@desy.de>, Katerina Lipka <katerina.lipka@desy.de>
// transition from pole to MSbar scheme by S. Moch (private communication)
// Modified on 2021-12-04 by T. Mäkelä (toni.makela@cern.ch): 
//   fixed accounting for LO and NLO in the MSBAR numerical stencil, made the
//   approach generic and generalizable to further schemes if need be. Included
//   s-channel and W+t final state processes in addition to previous t-channel.
  */

class IHathorGenericIntegrator;
template <class Integrand> class HathorGenericIntegrator;
class SgTop;
class HathorSgTopT;
class HathorSgTopS;
class HathorSgTopWt;
class HathorPdfxFitter;
class ReactionHathorSingleTop : public ReactionTheory
{
public:
    ReactionHathorSingleTop();
    ~ReactionHathorSingleTop();
    
    //vector<double> asFactors(SgTop *XS, double muOLD, double muNEW);
    vector<double> asFactors(int msMass, double nfl, int orderI, IHathorGenericIntegrator *XS, double muOLD, double muNEW);
    
    virtual string getReactionName() const { return  "HathorSingleTop" ;};
    virtual void initTerm(TermData *td) override final;
    virtual void atStart();
    virtual void compute(TermData *td, valarray<double> &val, map<string, valarray<double> > &err);
    virtual void atIteration();
protected:
    virtual int parseOptions(){ return 0;};
  
    HathorGenericIntegrator<HathorSgTopT>* _hathorT;
    HathorGenericIntegrator<HathorSgTopS>* _hathorS;
    HathorGenericIntegrator<HathorSgTopWt>* _hathorWt;
    
    HathorPdfxFitter* _pdf;
    int* _rndStore;
  
    std::string _integrator;
    const double* _sin2thW;
    const double* _alphaem;
    double _ckm[3][3];
    //Flags for which processes to include in the computation
    std::map<int, int> _tchannel;
    std::map<int, int> _schannel;
    std::map<int, int> _Wtchannel;
    std::map<int, int> _antiquark;
    std::map<int, const double* > _mtopPerInstance;
    std::map<int, double > _mrPerInstance;
    std::map<int, double > _mfPerInstance;
    std::map<int, std::string > _orderPerInstance;
    std::map<int, int > _MsbarPerInstance;
    std::map<int, int > _precisionLevelPerInstance;
    std::map<int, double > _sqrtSPerInstance;
    std::map<int, int > _ppbarPerInstance;

    std::map<std::string, std::pair<double, std::vector<TermData*> > > _convolved;
    std::vector<std::string> _convolved_vector_of_keys; // for parallel execuion

    // constants
    static constexpr double pi = 3.141592653589793;
    static constexpr double z2 = 1.644934066848226;
    static constexpr double z3 = 1.202056903159594;
    static constexpr double ln2= 0.693147180559945;
};
