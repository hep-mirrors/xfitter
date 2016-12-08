 
/*
   @file ReactionfastNLO.cc
   @date 2016-12-06
   @author  AddReaction.py
   Created by  AddReaction.py on 2016-12-06
*/

#include "ReactionfastNLO.h"
//#include "FastNLOxfReactionTheory.h" // todo
#include "fastnlotk/fastNLOReader.h"


// the class factories
extern "C" ReactionfastNLO* create() {
  return new ReactionfastNLO();
}


// Initialize at the start of the computation
void ReactionfastNLO::initAtStart(const string &s)
{
   // for testing purposes: create an object here
   fastNLOTable* fnlo = new fastNLOTable();
   fnlo->PrintHeader();
}

// Main function to compute results at an iteration
int ReactionfastNLO::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
   return 0;
}

