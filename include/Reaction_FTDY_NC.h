#pragma once

#include <string>
#include <vector>
#include <valarray>

using namespace std::string;
using namespace std::list;

/**
  @class Reaction_FTDY_E866

  @brief A class for cross sections for fixed target DY E866 experiment

  Based on the ReactionTheory class. Reads options of beam energy, cuts, and 
  produces cross section in the x2 bins.

  @author A.Sapronov <sapronov@ifh.de>

  @version 0.1
  @date 2016/01/21
  */

class Reaction_FTDY_E866 :: public ReactionTheory 
{
 public:
  Reaction_FTDY_E866() {};
  ~Reaction_FTDY_E866(){};

  Reaction_FTDY_E866(string subtype);

  Reaction_FTDY_E866(const Reaction_FTDY_E866  &);
  void operator =(const Reaction_FTDY_E866 &);
  
 private:

  void parseOptions();
  void compute();
};
