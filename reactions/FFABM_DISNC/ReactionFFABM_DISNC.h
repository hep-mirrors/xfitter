
#pragma once

#include "ReactionBaseDISNC.h"

/**
  @class' ReactionFFABM_DISNC

  @brief A wrapper class for FFABM_DISNC reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-09-29
  */

class ReactionFFABM_DISNC : public ReactionBaseDISNC
{
private:
  typedef ReactionBaseDISNC Super;

public:
  ReactionFFABM_DISNC(){};

public:
  virtual string getReactionName() const override { return "FFABM_DISNC"; };
  void virtual atStart() override final;
  virtual void initTerm(TermData *td) override final;
  virtual void atIteration() override final;

protected:
  virtual void F2 BASE_PARS override final;
  virtual void FL BASE_PARS override final;
  virtual void xF3 BASE_PARS override final;

private:
  map<int, valarray<double>> _f2abm;
  map<int, valarray<double>> _flabm;
  map<int, valarray<double>> _f3abm;

  // parameters initialised at iteration
  // (pointers for those parameters which can change at each iteration)
  const double* _mcPtr;
  const double* _mbPtr;
  const double* _mzPtr;
  const double* _sin2thwPtr;

  void calcF2FL(unsigned dataSetID);

  double apply_tmc(double& f2, double& fl, double& f3, const int flag_flavour, const std::valarray<double>& q2, const std::valarray<double>& x,
    const int ncflag, const int charge, const double polarity, const double cos2thw, const size_t i);
  struct integration_params {
    std::valarray<double> q2;
    int i;
    int ncflag;
    int charge;
    double polarity;
    double cos2thw;
    const double* _sin2thwPtr;
    const double* _mzPtr;
    int flag_calc_fl;
    int flag_flavour;
    int order;
  };
  int _flag_tmc;
};
