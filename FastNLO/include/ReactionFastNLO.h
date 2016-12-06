#pragma once

#include <string>
#include <vector>
#include <valarray>

#include "ReactionTheory.h"

using std::string;
using std::list;

/**
  @class ReactionFastNLO

  @brief A wrapper class for FastNLO 

  Based on the ReactionTheory class. 

  @version 0.2
  @date 2016/12/17
  */

class ReactionFastNLO : public ReactionTheory 
{
 public:
  ReactionFastNLO() {};
  ~ReactionFastNLO(){};

  ReactionFastNLO(const ReactionFastNLO  &){};
  ReactionFastNLO & operator =(const ReactionFastNLO &r){return *(new ReactionFastNLO(r));}; // this causes a memory leak!

  string getReactionName() const {return string("fastNLO");};
  
 public:
  int compute();
 private:

  void initAtStart(const string &); 
  int parseOptions(){};
};
