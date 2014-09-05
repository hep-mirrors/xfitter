// Author: Daniel Britzger
// DESY, 29/01/2014

/**
   fastNLOInterpolOneNode

   No interpolation is performed, but all
   values are stored at one single node.
   The 'x'-value of the node is calculated
   from the H-function.
*/

#ifndef __fastNLOInterpolOneNode__
#define __fastNLOInterpolOneNode__

#include "speaker.h"
#include <string>
#include <vector>
#include <utility>
#include "fastNLOInterpolBase.h"

using namespace std;

class fastNLOInterpolOneNode :  public fastNLOInterpolBase {
   
public:

   fastNLOInterpolOneNode(double min, double max);
   ~fastNLOInterpolOneNode(void);
   
   //   vector<pair<int,double> > CalcNodeValues(double val);
   void CalcNodeValues(vector<pair<int,double> >& nodes, double val);

protected:


private:
   vector<pair<int,double> > fDummyNode;

};


#endif
