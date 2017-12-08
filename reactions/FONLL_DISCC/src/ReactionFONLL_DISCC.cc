 
/*
   @file ReactionFONLL_DISCC.cc
   @date 2017-11-29
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-11-29
*/

#include "ReactionFONLL_DISCC.h"

// the class factories
extern "C" ReactionFONLL_DISCC* create() {
  return new ReactionFONLL_DISCC();
}


// Initialize at the start of the computation
int ReactionFONLL_DISCC::initAtStart(const string &s)
{
  return 0;
}

// Main function to compute results at an iteration
int ReactionFONLL_DISCC::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
  return 0;
}

