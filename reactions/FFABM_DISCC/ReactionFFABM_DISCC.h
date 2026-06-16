
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

struct ReactionFFABM_DISCC;
struct nomad_integration_params {
  int intvar;
  double val;
  unsigned dataSetID;
  const BaseDISCC::ReactionData* rd;
  ReactionFFABM_DISCC* reaction;
  const double* br0;
  const double* br1;
  double mnucl;
  int nomad_scaleq2mw2;
  int nomad_scalesemilepbr;
};

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
  std::map<unsigned, int> _order; // term dependent, allows adjustment for specific data sets
  std::map<unsigned, bool> _msbarmin;
  std::map<unsigned, int> _ordfl;
  std::map<unsigned, int> _orderHQ;
  std::vector<const double*> _ckm; // CKM matrix
  bool _need_pdffillgrid;

  void calcF2FL(int dataSetID);
  void calc_point(const double q2, const double x, const int dataSetID, const BaseDISCC::ReactionData *rd, double& f2out, double& flout, double& f3out);
  //double calc_point_strfun(const BaseDISCC::ReactionData* rd, const BaseDISCC::dataType ftype, const BaseDISCC::dataFlav flav, const double q2, const double x, const int dataSetID, const int order, const int charge, const double f2c=0.);
  void calc_integral_cuba(const int intvar, const double val, const int dataSetID, const BaseDISCC::ReactionData *rd, double& xsec_out, const int nomad_scaleq2mw2, const int nomad_scalesemilepbr, const double nomad_epsrel, const int nomad_verbose);
  static int integrate_nomad_cubareal(const int* ndim, const cubareal* x, const int *ncomp, cubareal* ff, void *userdata);
  static int integrate_nomad(const int* ndim, const double* x, const int *ncomp, double* ff, void *userdata);
  static long unsigned _ncalls_nomad;
  static double integrand_sd2(double x[]);
  static void integrand_sd2_region(int ll, double* xx, double* aa, double* bb);
  int _sd2_nomad_var;
  nomad_integration_params make_integration_params(const int intvar, const double val, const int dataSetID, const BaseDISCC::ReactionData *rd, const int nomad_scaleq2mw2, const int nomad_scalesemilepbr);
  nomad_integration_params _sd2_nomad_pars;
  static nomad_integration_params _sd2_nomad_pars_static;
  map<int, int> _ncpu;
  double combine_flavours(const BaseDISCC::ReactionData* rd, const double f, const double fc, const double fb);
  double _q2mincomp;
};

