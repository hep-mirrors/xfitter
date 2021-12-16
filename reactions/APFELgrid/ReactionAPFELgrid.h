#pragma once
#include"ReactionTheory.h"
/**
  @class' ReactionAPFELgrid

  @brief A wrapper class for APFELgrid reaction
  */
class ReactionAPFELgrid:public ReactionTheory{
public:
  virtual string getReactionName()const{return"APFELgrid";};
  virtual void initTerm(TermData*)override final;
  virtual void freeTerm(TermData*)override final;
  virtual void compute(TermData*,valarray<double>&,map<string,valarray<double> >&)override final;
};
