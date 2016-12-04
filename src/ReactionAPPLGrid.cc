/*!
 @file CommonGrid.cc
 @date Tue March 25 2014
 @author Andrey Sapronov <Andrey.Sapronov@cern.ch>

 Contains CommonGrid class member function implementations.
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>
#include <float.h>

#include <string.h>

#include "ReactionAPPLGrid.h"

using std::string;

// the class factories

extern "C" ReactionAPPLGrid* create() {
  return new ReactionAPPLGrid();
}


int ReactionAPPLGrid::compute()
{
}

void ReactionAPPLGrid::initAtStart(const string &s)
{
  std::cout << s << std::endl;
}
