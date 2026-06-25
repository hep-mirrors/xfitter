
#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionHathor

  @brief A wrapper class for Hathor reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-08-07
  */

//class Hathor;
template <class Integrand> class HathorGenericIntegrator;
class HathorLikeSgTop;
class HathorPdfxFitter;

class ReactionHathor : public ReactionTheory
{
public:
  ReactionHathor();

  ~ReactionHathor();

  //    ~ReactionHathor(){};
  //    ~ReactionHathor(const ReactionHathor &){};
  //    ReactionHathor & operator =(const ReactionAHathor &r){return *(new ReactionHathor(r));};

public:
  virtual string getReactionName() const { return  "Hathor" ;};
  virtual void initTerm(TermData *td) override final;
  virtual void atStart();
  virtual void compute(TermData *td, valarray<double> &val, map<string, valarray<double> > &err);
  virtual void atIteration();
protected:
  HathorGenericIntegrator<HathorLikeSgTop>* _hathor;
  HathorPdfxFitter* _pdf;
  int* _rndStore;
  bool _init = false;
  std::string _integrator;
  std::map<int, const double* > _mtopPerInstance;
  std::map<int, double > _mrPerInstance;
  std::map<int, double > _mfPerInstance;
  std::map<int, std::string > _orderPerInstance;
  std::map<int, int > _NfPerInstance;
  std::map<int, int > _MsbarPerInstance;
  std::map<int, int > _precisionLevelPerInstance;
  std::map<int, double > _sqrtSPerInstance;
  std::map<int, int > _ppbarPerInstance;
  std::map<std::string, std::pair<double, std::vector<TermData*> > > _convolved;
  std::vector<std::string> _convolved_vector_of_keys; // for parallel execuion
};

