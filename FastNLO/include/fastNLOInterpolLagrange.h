// Author: Daniel Britzger
// DESY, 29/01/2014

/**
   fastNLOInterpolLagrange

   Interpolation routines for lagrange interpolation
   of second order polynomials.
*/

#ifndef __fastNLOInterpolLagrange__
#define __fastNLOInterpolLagrange__

#include "speaker.h"
#include <string>
#include <vector>
#include <utility>
#include "fastNLOInterpolBase.h"

using namespace std;

class fastNLOInterpolLagrange : public fastNLOInterpolBase {
   
public:

   fastNLOInterpolLagrange(double min, double max);
   ~fastNLOInterpolLagrange(void);
   
   //   vector<pair<int,double> > CalcNodeValues(double val);
   void CalcNodeValues(vector<pair<int,double> >& nodes, double val);

protected:


private:


};


#endif
