
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
  double _hqscale1in;
  double _hqscale2in;
  bool _msbarmin;
  int _ordfl;
  std::map<unsigned, int> _order; // term dependent
  int _kschemepdfin;
  int _FLAG_FAST = 0;
  std::map<unsigned, int> _orderHQ; // allows adjusting HQ order for specific data sets

  void calcF2FL(int dataSetID);
  void calc_point(const double q2, const double x, const int dataSetID, const BaseDISCC::ReactionData *rd, double& f2out, double& flout, double& f3out);
  //double calc_point_strfun(const BaseDISCC::ReactionData* rd, const BaseDISCC::dataType ftype, const BaseDISCC::dataFlav flav, const double q2, const double x, const int dataSetID, const int order, const int charge, const double f2c=0.);
  double calc_point_strfun(const BaseDISCC::ReactionData* rd, const BaseDISCC::dataType ftype, const BaseDISCC::dataFlav flav, const double q2, const double x, const int dataSetID, const int order, const int charge, const double* f2c=nullptr);
  void calc_integral(const int intvar, const double val, const int dataSetID, const BaseDISCC::ReactionData *rd, double& xsec_out);
  static int integrate_nomad(const int* ndim, const cubareal* x, const int *ncomp, cubareal* ff, void *userdata);
  map<int, int> _ncpu;
  double combine_flavours(const BaseDISCC::ReactionData* rd, const double f, const double fc, const double fb);
};

