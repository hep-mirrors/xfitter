// Author: Daniel Britzger
// DESY, 29/01/2014

/**
   fastNLOInterpolLinear

   Interpolation routines for linear interpolation.
*/

#ifndef __fastNLOInterpolLinear__
#define __fastNLOInterpolLinear__

#include "speaker.h"
#include <string>
#include <vector>
#include <utility>
#include "fastNLOInterpolBase.h"

using namespace std;

class fastNLOInterpolLinear : public fastNLOInterpolBase {
   
public:

   fastNLOInterpolLinear(double min, double max);
   ~fastNLOInterpolLinear(void);
   
   //   vector<pair<int,double> > CalcNodeValues(double val);
   void CalcNodeValues(vector<pair<int,double> >& nodes, double val);

protected:


private:


};


#endif
