
#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionKMatrix

  @brief A wrapper class for KMatrix reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2018-08-03
  */

class ReactionKMatrix:public ReactionTheory{
public:
  ReactionKMatrix(){};
  virtual string getReactionName() const override {return"KMatrix";}
  virtual void initTerm(TermData*) override final;
  virtual void freeTerm(TermData*) override final;
  virtual void compute(TermData*,valarray<double>&val,map<string,valarray<double> >&errors) override final;
};

