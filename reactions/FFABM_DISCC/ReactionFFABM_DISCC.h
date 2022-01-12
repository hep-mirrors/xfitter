
#pragma once

#include "ReactionBaseDISCC.h"

/**
  @class' ReactionFFABM_DISCC

  @brief A wrapper class for FFABM_DISCC reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-10-09
  */

class ReactionFFABM_DISCC : public ReactionBaseDISCC
{
private:
  typedef ReactionBaseDISCC Super;
public:
  ReactionFFABM_DISCC(){};
  virtual string getReactionName() const override { return  "FFABM_DISCC" ;};
  void virtual atStart() override final;
  virtual void initTerm(TermData *td) override final;
  virtual void atIteration() override final;

protected:
  virtual valarray<double> F2(TermData *td) override final;
  virtual valarray<double> FL(TermData *td) override final;
  virtual valarray<double> xF3(TermData *td) override final;

private:
  map <int,valarray<double> > _f2abm;
  map <int,valarray<double> > _flabm;
  map <int,valarray<double> > _f3abm;

  // parameters initialised at iteration
  // (pointers for those parameters which can change at each iteration)
  const double* _mcPtr;
  const double* _mbPtr;
  const double* _mzPtr;
  const double* _sin2thwPtr;

  void calcF2FL(int dataSetID);
};

