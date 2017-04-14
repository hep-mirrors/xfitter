 
/*
   @file ReactionTensorPomeron.cc
   @date 2017-04-14
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-04-14
*/

#include "ReactionTensorPomeron.h"

// the class factories
extern "C" ReactionTensorPomeron* create() {
  return new ReactionTensorPomeron();
}


// Initialize at the start of the computation
int ReactionTensorPomeron::initAtStart(const string &s)
{
  return 0;
}

// Main function to compute results at an iteration
int ReactionTensorPomeron::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
  return 0;
}

