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
  virtual void atStart() override final;
  virtual void atIteration() override final;
  virtual void initTerm(TermData*)override final;
 
 protected:
  virtual valarray<double> F2(TermData*td) override;
  virtual valarray<double> FL(TermData*td) override;
  virtual valarray<double> xF3(TermData*td) override;

 private:
  map <unsigned,valarray<double>> _f2fonll;
  map <unsigned,valarray<double>> _flfonll;
  map <unsigned,valarray<double>> _f3fonll;
  map <unsigned,TermData*> _dsIDs;
  //! Allow for non-apfelff evolution:
  bool _non_apfel_evol{false};  
};

