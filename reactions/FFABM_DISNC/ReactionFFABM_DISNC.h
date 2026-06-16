
#pragma once

#include "ReactionBaseDISNC.h"
#include "ForkPool.h"
//#include "cuba.h"

/**
  @class' ReactionFFABM_DISNC

  @brief A wrapper class for FFABM_DISNC reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-09-29
  */

class ReactionFFABM_DISNC : public ReactionBaseDISNC
{
  friend DIS_TMC; // TODO modify DIS_TMC to call public methods only
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
  double _hqscale1in;
  double _hqscale2in;
  std::map<unsigned, int> _order; // term dependent, allows adjustment for specific data sets
  std::map<unsigned, bool> _msbarmin;
  std::map<unsigned, int> _ordfl;
  std::map<unsigned, int> _orderHQ;
  bool _need_pdffillgrid;

  void calcF2FL(unsigned dataSetID);

  double apply_tmc(const int method, double& f2, double& fl, double& f3, const int flag_flavour, const std::valarray<double>& q2, const std::valarray<double>& x,
    const int ncflag, const int charge, const double polarity, const double cos2thw, const size_t i);
  //static int Integrand_Cuhre(const int* ndim, const cubareal* x, const int *ncomp, cubareal* ff, void *userdata);
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
  };
  map<int, int> _flag_tmc;
  map<int, int> _tmc_integration_method;
  map<int, double> _tmc_xmin;
  map<int, double> _tmc_logxlogq2min;
  const double* _tmc_mpr;
  map<int, int> _ncpu;
  static double combine_flavours(const ReactionBaseDISNC::dataFlav flav, const double f, const double fc, const double fb);
  double _q2mincomp;
  map<int, bool> _calcf2fldone;

  struct DataPoint {
    int datasetID;
    int i;
    ReactionBaseDISNC::dataFlav flav;
    int ord;
    int ordHQ;
    int ordFL;
    int msbarmin;
    double charge;
    double polar;
    const double* sin2thetaWPtr;
    double mz;
    double q2;
    double x;
    DIS_NUKE* nuke;
    DIS_TMC* tmc;
    DIS_HT* ht;
    double f2;
    double fl;
    double f3;
    void calc();
  };
  ForkPool::TaskDistribution _task_distr; // 0 is chunky, 1 is cyclic
  // todo optimize memory
  std::map<int, std::vector<DataPoint> > _data_points;
  std::map<std::string, std::vector<DataPoint> > _grouped_data_points;
};