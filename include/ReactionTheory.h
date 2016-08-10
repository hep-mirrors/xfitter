#pragma once

#include <string>
#include <vector>
#include <valarray>

using std::string;
using std::vector;
using std::valarray;

/**
  @class ReactionTheory

  @brief A base class manages for reaction theories

  It provides an interface wich must present in the derived classes

  @author A.Sapronov <sapronov@ifh.de>

  @version 0.1
  @date 2016/01/21
  */

class ReactionTheory 
{
 public:
  ReactionTheory() {};
  ~ReactionTheory();

  ReactionTheory(string subtype) : _subtype(subtype) {};

  ReactionTheory(const ReactionTheory &);
  ReactionTheory & operator =(const ReactionTheory &);

 public:
  void setOptions(const string &reaction_options) { _ro = reaction_options; this->parseOptions(); };
  void setBinning(vector<int> &binFlags, vector<vector<double> > &dsBins){ _binFlags=binFlags; _dsBins = dsBins; } ;
  void resultAt(valarray<double> *val){ _val = val; };
  
  virtual int compute() = 0;
 protected:

  virtual int parseOptions() = 0;

 protected:
  string _subtype;
  valarray<double> *_val;
  string _ro;
  vector<int> _binFlags;
  /// dataset bins
  vector<vector<double> > _dsBins;
};
