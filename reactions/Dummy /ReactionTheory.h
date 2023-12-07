#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionDummy

  @brief A wrapper class for Dummy reaction

  Based on the ReactionTheory class. 

  @version 1.0
  @date 2023-12-15
  */

class ReactionDummy:public ReactionTheory{
public:
  virtual string getReactionName() const override {return "Dummy";};
  virtual void compute(TermData*,valarray<double>&,map<string,valarray<double> >&) override final;
};
