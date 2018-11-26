#pragma once

#include "ReactionBaseDISCC.h"

/*
 *   @class' ReactionFONLL_DISCC
 *
 *  @brief A wrapper class for FONLL_DISCC reaction
 *
 *  Based on the ReactionTheory class. Reads options produces 3d cross section.
 *
 *  @version 0.1
 *  @date 2017-11-29
 */

class ReactionFONLL_DISCC : public ReactionBaseDISCC
{
 public:
  ReactionFONLL_DISCC() {};
  //~ReactionFONLL_DISCC() {};
  //~ReactionFONLL_DISCC(const ReactionFONLL_DISCC &) {};
  //ReactionFONLL_DISCC & operator = (const ReactionAFONLL_DISCC &r) { return *(new ReactionFONLL_DISCC(r)); };

  virtual string getReactionName() const { return "FONLL_DISCC"; };
  int atStart(const string &);
  virtual void initAtIteration() override;

 protected:
  virtual void F2  BASE_PARS override;
  virtual void FL  BASE_PARS override;
  virtual void xF3 BASE_PARS override;

 private:
  map <int,valarray<double>> _f2fonll;
  map <int,valarray<double>> _flfonll;
  map <int,valarray<double>> _f3fonll;
};

