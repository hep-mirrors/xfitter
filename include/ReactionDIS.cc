#pragma once

#include <string>
#include <vector>
#include <valarray>

#include "ReactionTheory.h"

using namespace std::string;
using namespace std::list;

/**
  @class ReactionDIS

  @brief A class for cross sections for HERA DIS NC experiment

  Based on the ReactionTheory class. Reads options of beam energy, cuts, and 
  produces cross section in the x2 bins.

  @author A.Sapronov <sapronov@ifh.de>

  @version 0.1
  @date 2016/05/17
  */

class ReactionDIS :: public ReactionTheory 
{
 public:
  ReactionDIS() {};
  ~ReactionDIS(){};

  ReactionDIS(string subtype);

  ReactionDIS(const ReactionDIS  &);
  void operator =(const ReactionDIS &);
  
 private:

  void parseOptions();
  void compute();
};
