
#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionFractal_DISNC

  @brief A wrapper class for Fractal_DISNC reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-05-15
  */

class ReactionFractal_DISNC:public ReactionTheory{
public:
  virtual string getReactionName() const override {return"Fractal_DISNC";};
  virtual void compute(TermData*,valarray<double>&,map<string,valarray<double> >&) override final;
};

