
#pragma once

#include "ReactionBaseDISNC.h"
#include "ReactionBaseDISCC.h"
#include "ForkPool.h"
//#include "cuba.h"

class ReactionFFABM_DISNC_CC : public ReactionBaseDISCC
//class ReactionFFABM_DISNC_CC : public ReactionBaseDISNC, public ReactionBaseDISCC
{
private:
  typedef ReactionBaseDISCC Super;

public:
  ReactionFFABM_DISNC_CC(){};

public:
  virtual string getReactionName() const override { return "FFABM_DISNC_CC"; };
  void virtual atStart() override final;
  virtual void initTerm(TermData *td) override final;
  virtual void atIteration() override final;

protected:
  virtual valarray<double> F2(TermData *td) override final;
  virtual valarray<double> FL(TermData *td) override final;
  virtual valarray<double> xF3(TermData *td) override final;

private:
  map<int, valarray<double>> _f2abm;
  map<int, valarray<double>> _flabm;
  map<int, valarray<double>> _f3abm;

  // pointers to those parameters which can change at each iteration
  const double* _mcPtr;
  const double* _mbPtr;

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
  static double combine_flavours(const ReactionBaseDISCC::dataFlav flav, const double f, const double fc, const double fb);
  double _q2mincomp;

  struct DataPoint {
    TermData* td;
    int datasetID;
    int i;
    ReactionBaseDISCC::dataFlav flav;
    int ord;
    int ordHQ;
    int ordFL;
    int msbarmin;
    double charge;
    double polar;
    const double* sin2thetaWPtr;
    const double* mz;
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
  std::map<std::string, std::vector<DataPoint> > _grouped_data_points;
};
