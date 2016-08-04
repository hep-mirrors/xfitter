#pragma once

#include <string>
#include <vector>
#include <valarray>

#include "ReactionTheory.h"

using namespace std::string;
using namespace std::list;

/**
  @class ReactionDIS

  @brief A wrapper class for cross DIS NC sections of HERA

  Based on the ReactionTheory class. Reads options produces 3d cross section.

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
