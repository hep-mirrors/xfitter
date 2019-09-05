
#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionUbarDbar

  @brief A wrapper class for UbarDbar reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2019-09-05
  */

class ReactionUbarDbar : public ReactionTheory
{
  public:
    ReactionUbarDbar(){};


  public:
    virtual string getReactionName() const { return  "UbarDbar" ;};
    virtual void compute(TermData* td, valarray<double> &val, map<string, valarray<double> > &err) override final;
};

