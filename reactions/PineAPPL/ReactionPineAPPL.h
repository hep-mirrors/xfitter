#pragma once
#include"ReactionTheory.h"
/**
  @class' ReactionPineAPPL

  @brief A wrapper class for PineAPPL reaction
  */
class ReactionPineAPPL:public ReactionTheory{
  std::map<std::string, vector<double> > _convolved;
public:
  virtual string getReactionName() const override {return"PineAPPL";};
  virtual void initTerm(TermData*) override final;
  virtual void freeTerm(TermData*) override final;
  virtual void compute(TermData*,valarray<double>&,map<string,valarray<double> >&) override final;
  virtual void atIteration();
};
//NOTE: parameter "energy" is only read at start, once, and therefore cannot be fitted
