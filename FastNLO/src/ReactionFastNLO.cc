/*!
 @file ReactionFastNLO.cc
 @date Dec 2016

 Contains implementation of ReactionFastNLO based on ReactionTheory
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>
#include <float.h>
#include <string.h>

#include <ReactionFastNLO.h>

using std::string;

// the class factories

extern "C" ReactionFastNLO* create() {
  return new ReactionFastNLO();
}


int ReactionFastNLO::compute()
{
}

void ReactionFastNLO::initAtStart(const string &s)
{
  std::cout << s << std::endl;
}
