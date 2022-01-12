#pragma once
#include"ReactionTheory.h"
/**
  @class' ReactionAPPLgrid

  @brief A wrapper class for APPLgrid reaction
  */
class ReactionAPPLgrid:public ReactionTheory{
public:
  virtual string getReactionName() const override {return"APPLgrid";};
  virtual void initTerm(TermData*) override final;
  virtual void freeTerm(TermData*) override final;
  virtual void compute(TermData*,valarray<double>&,map<string,valarray<double> >&) override final;
};
//NOTE: parameter "energy" is only read at start, once, and therefore cannot be fitted
