
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
  virtual string getReactionName() const { return  "TensorPomeron" ;};
  virtual int initAtStart(const string &s) override;
  virtual void initAtIteration() override; 
  virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err) override;
  virtual void setDatasetParamters( int dataSetID, map<string,string> parsReaction,  map<string,double> parsDataset) override;
 
 private:
  map <int, int> _npoints;
  double _e0, _e1;  ///< pomeron slopes
  double _mp, _alpha_em;
  double _alphaP, _beta;
  double _m02, _m12;

  // the main fit target:
  double _epsilon0, _epsilon1;

  // store also kin vars:
  map <int, valarray<double> > _x,_q2,_y, _W2,_delta,_Yp, _f, _pq;
  
  // and derrived too:
  map <int, valarray<double> > _b0,_b1,_qa0,_qa1;

  // spline functions
  vector<double> _s0bn, _s1bn, _s0rn, _s1rn; ///< spline knots
  tk::spline _s0b, _s1b, _s0r, _s1r;

  void sigma_L (int dataSetID, valarray<double>& sL);
  void sigma_LT(int dataSetID, valarray<double>& sLT);

  // Q2 dependences
  const double b0q2(double q2){ return exp( _s0b(q2+_m02) ); }
  const double b1q2(double q2){ return exp( _s1b(q2+_m12) ); }

  const double r0q2(double q2){ return exp( _s0r(q2) ); }
  const double r1q2(double q2){ return exp( _s1r(q2) ); }

  const double q2a0(double q2) { return b0q2(q2)/(1. + 1./r0q2(q2)); }
  const double q2a1(double q2) { return b1q2(q2)/(1. + 1./r1q2(q2)); }
  
};

