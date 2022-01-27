
#pragma once

#include "ReactionTheory.h"
#include "spline.h"

/**
  @class' ReactionTensorPomeron

  @brief A wrapper class for TensorPomeron reaction 

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-04-14
  */

class ReactionTensorPomeron : public ReactionTheory
{
 public:
  ReactionTensorPomeron(){};
 public:
  virtual string getReactionName() const override { return  "TensorPomeron" ;};
  virtual void atIteration() override; 
  virtual void compute(TermData* td, valarray<double> &val, map<string, valarray<double> > &err) override;
  virtual void initTerm(TermData* td) override;
  virtual void atFCN3() override;
  virtual void atMakeErrorBands (int ivector) override; 
  const valarray<double> *GetBinValues(TermData *td, const string &binName);

 private:
  vector<unsigned> _dsIDs;
  map <int, int> _npoints;
  double _mp, _alpha_em;
  double _alphaP, _beta;
  double _m02, _m12;
  double _convFact;

  // hard and soft pomeron slopes:
  double _epsilon0, _epsilon1;

  map <int,bool> _isReduced;

  /// Flag for spline param. for R
  int _splineR;
  double _mr, _r0, _r1, _delta0, _delta1;

  /// Reggeon
  double _alphaR, _betaR, _epsilonR, _alphaIR, _c0ir, _c1ir;

  // store also kin vars:
  map <int, valarray<double> > _x,_q2,_y, _W2,_delta,_Yp, _f, _pq;
  
  // and derrived too:
  map <int, valarray<double> > _b0,_b1,_qa0,_qa1, _b2, _qa2; 

  // spline functions
  vector<double> _s0bn, _s1bn, _s0rn, _s1rn; ///< spline knots
  tk::spline _s0b, _s1b, _s0r, _s1r;

  void sigma_L (int dataSetID, valarray<double>& sL);
  void sigma_LT(int dataSetID, valarray<double>& sLT);

  // Q2 dependences
  const double b0q2(double q2){ return exp( _s0b(log(q2+_m02)) ); }
  const double b1q2(double q2){ return exp( _s1b(log(q2+_m12)) ); }

  // Reggeon
  const double b2q2(double q2){ return exp(_c0ir - q2/_c1ir);  }

  const double r0q2(double q2);
  const double r1q2(double q2);
  const double r2q2(double q2) {return 0.; } //!< Reggeon FL=0.

  const double q2a0(double q2) { return 0.5*b0q2(q2)/(1. + 1./r0q2(q2)); }
  const double q2a1(double q2) { return 0.5*b1q2(q2)/(1. + 1./r1q2(q2)); }
  const double q2a2(double q2) { return 0.5*b2q2(q2)/(1. + 1./r2q2(q2)); }

  void writeOut(const std::string& file);   ///< Store b and a functions.
  vector<double> getSplinePar(const std::string& vn); ///< Helper to decode variable pars.
  
};

