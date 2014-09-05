// Author: Daniel Britzger
// DESY, 28/06/2013

#ifndef __fastNLOInterpolCatmullRom__
#define __fastNLOInterpolCatmullRom__

#include "speaker.h"
#include <string>
#include <vector>
#include <utility>
#include "fastNLOInterpolBase.h"

using namespace std;

class fastNLOInterpolCatmullRom : public fastNLOInterpolBase {

public:

   fastNLOInterpolCatmullRom(double min, double max);
   ~fastNLOInterpolCatmullRom(void);

   //   vector<pair<int,double> > CalcNodeValues(double val);
   void CalcNodeValues(vector<pair<int,double> >& nodes, double val);

protected:


private:


};
#endif
