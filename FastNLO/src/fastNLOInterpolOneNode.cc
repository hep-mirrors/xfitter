// Author: Daniel Britzger
// DESY, 29/01/2014

#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include "fastNLOInterpolOneNode.h"

using namespace std;

//______________________________________________________________________________
fastNLOInterpolOneNode::fastNLOInterpolOneNode(double min, double max) : fastNLOInterpolBase(min,max,1) {
   debug["fastNLOInterpolOneNode"]<<"New fastNLOInterpolOneNode instance."<<endl;
   fDummyNode.resize(1);
   fDummyNode[0] = make_pair(0,1);
   if (fLastGridPointWasRemoved ) 
      warn["fastNLOInterpolOneNode"]<<"Last grid point cannot be removed, since there is only one point."<<endl;
}


//______________________________________________________________________________
fastNLOInterpolOneNode::~fastNLOInterpolOneNode(void) {
}


//______________________________________________________________________________
void fastNLOInterpolOneNode::CalcNodeValues(vector<pair<int,double> >& nodes, double x) {
   nodes=fDummyNode;
}

