#pragma once
#include "ReactionTheory.h"

/**
  @class' ReactionBaseDISCC

  @brief A wrapper class for BaseDISCC reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-10-05
  */
class ReactionBaseDISCC : public ReactionTheory
{
  public:
    ReactionBaseDISCC(){};
  public:
    virtual string getReactionName()const{return"BaseDISCC";};
    virtual void atStart()override final;
    virtual void initTerm(TermData*)override final;
    virtual void compute(TermData*,valarray<double>&val,map<string,valarray<double> >&errors)override final;
  protected:
    double _Gf;
    double _convfac;
};

