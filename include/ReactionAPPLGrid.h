#pragma once

#include <string>
#include <vector>
#include <valarray>

#include "ReactionTheory.h"

using std::string;
using std::list;

/**
  @class ReactionAPPLGrid

  @brief A wrapper class for APPLGrid 

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2016/05/17
  */

class ReactionAPPLGrid : public ReactionTheory 
{
 public:
  ReactionAPPLGrid() {};
  ~ReactionAPPLGrid(){};

  ReactionAPPLGrid(const ReactionAPPLGrid  &){};
  ReactionAPPLGrid & operator =(const ReactionAPPLGrid &r){return *(new ReactionAPPLGrid(r));};
  
 public:
  int compute();
 private:

  void initAtStart(const string &); 
  int parseOptions(){};
};
