#pragma once

#include <string>
#include <vector>
#include <valarray>

#include "ReactionTheory.h"

using std::string;
using std::list;

/**
  @class ReactionDIS

  @brief A wrapper class for cross DIS NC sections of HERA

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @author A.Sapronov <sapronov@ifh.de>

  @version 0.1
  @date 2016/05/17
  */

class ReactionDIS : public ReactionTheory 
{
 public:
  ReactionDIS() {};
  ~ReactionDIS();

  ReactionDIS(const string &subtype);

  ReactionDIS(const ReactionDIS  &);
  ReactionDIS & operator =(const ReactionDIS &);
  
 public:
  int compute();
 private:

  int parseOptions();
};
