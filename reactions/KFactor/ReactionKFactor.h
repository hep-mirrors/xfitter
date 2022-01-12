
#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionKFactor

  @brief A wrapper class for KFactor reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-10-28
  */

class ReactionKFactor : public ReactionTheory{
public:
  ReactionKFactor(){};
  virtual string getReactionName() const override {return"KFactor";}
  virtual void initTerm(TermData*) override final;
  virtual void freeTerm(TermData*) override final;
  virtual void compute(TermData*,valarray<double>&val,map<string,valarray<double> >&errors) override final;
};

