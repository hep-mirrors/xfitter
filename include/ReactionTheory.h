#pragma once

#include <map>
#include <string>
#include <vector>
#include <valarray>

using std::map;
using std::string;
using std::vector;
using std::valarray;

typedef double (*pxFx)(double*, double*);
typedef double (*pZeroParFunc)();
typedef double (*pOneParFunc)(double*);
typedef double (*pTwoParFunc)(double*, double*);

/**
  @class ReactionTheory

  @brief A base class manages for reaction theories

  It provides an interface wich must present in the derived classes

  @author A.Sapronov <sapronov@ifh.de>

  @version 0.1
  @date 2016/01/21
  */

//class Evolution;

class ReactionTheory 
{
 public:
  ReactionTheory() {};
  ~ReactionTheory() {};

  ReactionTheory(const ReactionTheory &);
  ReactionTheory & operator =(const ReactionTheory &);

 public:
  virtual string getReactionName() const =0;
  virtual void initAtStart(const string &) =0;
  virtual void setxFitterParameters(map<string,double> &xfitter_pars) {*_xfitter_pars = xfitter_pars; };
  virtual void setEvolFunctions(double (*palpha_S)(double *) , map<string, pxFx> &) { alpha_S = palpha_S; };
  virtual void setExtraFunctions(map<string, pZeroParFunc>, map<string, pOneParFunc>, map<string, pTwoParFunc>) { };
  virtual void initAtIteration() {};
  virtual void setBinning(map<string,vector<double> > dsBins){ *_dsBins = dsBins; } ;
//  virtual void resultAt(valarray<double> *val){ _val = val; };
  
  virtual int compute(valarray<double> &val, map<string, valarray<double> > &err) = 0;
 protected:

  virtual int parseOptions() { return 0;};
  double (*alpha_S)(double *);

 protected:
  string _subtype;
  valarray<double> *_val;
  string _ro;
  /// dataset bins
  /// must contain 'binFlag' key
  map<string, vector<double> > *_dsBins;
  map<string, double > *_xfitter_pars;
};


typedef ReactionTheory * create_t();
