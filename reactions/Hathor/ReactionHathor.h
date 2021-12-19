
#pragma once

#include "ReactionTheory.h"
#include "Hathor.h"

/**
  @class' ReactionHathor

  @brief A wrapper class for Hathor reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-08-07
  */

class Hathor;
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
protected:
  virtual int parseOptions(){ return 0;};

  // this is map of key = dataset, value = pointer to Hathor instances,
  // one instance per one dataset
  //std::map<int, Hathor*> _hathorArray;

  //Global Hathor instance for all terms
  Hathor* _hathor;

  HathorPdfxFitter* _pdf;
  int* _rndStore;
  //double _mtop;
  std::map<int, std::shared_ptr<double> > _mtopPerInstance;
  std::map<int, std::shared_ptr<double> > _mrPerInstance;
  std::map<int, std::shared_ptr<double> > _mfPerInstance;
  // store term data for later access (never used)
  //map<unsigned, TermData*> _tdDS;
};

