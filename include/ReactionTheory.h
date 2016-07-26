#pragma once

#include <string>
#include <vector>
#include <valarray>

using namespace std::string;
using namespace std::list;

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
  ~ReactionTheory(){};

  ReactionTheory(string subtype) : _subtype(subtype) {};

  ReactionTheory(const ReactionTheory &);
  void operator =(const ReactionTheory &);

 public:
  void setOptions(const string &reaction_options) { _ro = reaction_options; this->parseOptions(); };
  void setBinning(vector<int> &binFlags, vector<vector<double> > &dsBins){ _binFlags=binFlags; _dsBins = dsBins; } ;
  void resultAt(valarray<double> *val){ _val = val; };
  
 protected:

  virtual void parseOptions() = 0;
  virtual void compute() = 0;

 protected:
  string _subtype;
  valarray<double> *_val;
  string _ro;
  vector<int> _binFlags;
  /// dataset bins
  vector<vector<double> > _dsBins;
};
