
#pragma once

#include "ReactionBaseDISCC.h"
#include "cuba.h"

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
  double _cos2thw;

  void calcF2FL(int dataSetID);
  void calc_point(const double q2, const double x, const int datasetID, const ReactionData *rd, double& f2out, double& flout, double& f3out);
  void calc_integral(const int intvar, const double val, const int datasetID, const ReactionData *rd, double& f2out, double& flout, double& f3out);
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
    double xi;
    int nt;
  };
  map<int, int> _flag_tmc;
  map<int, int> _flag_tmc_c;
  map<int, int> _flag_tmc_b;
  map<int, int> _tmc_integration_method;
  map<int, double> _tmc_xmin;
  map<int, double> _tmc_logxlogq2min;
  const double* _tmc_mpr;
  map<int, int> _ncpu;
  void combine_flavours(const int dataSetID, const double f, const double fc, const double fb, double& fout);
  double apply_tmc(const int method, double& f2, double& fl, double& f3, const bool flag_fl, const bool flag_f3, const int flag_flavour, const double q2, const double x, const int ncflag, const int charge, const double polarity, const double cos2thw, const int nt);
  static int Integrand_Cuhre(const int* ndim, const cubareal* x, const int *ncomp, cubareal* ff, void *userdata);
};

