#ifndef __fastNLOCoeffData__
#define __fastNLOCoeffData__

#include "fastNLOCoeffBase.h"
#include "fastNLOConstants.h"

using namespace std;

class fastNLOCoeffData : public fastNLOCoeffBase {

   friend class fastNLOTable;

public:
   fastNLOCoeffData(int NObsBin);
   fastNLOCoeffData(const fastNLOCoeffBase&);
   virtual ~fastNLOCoeffData(){;};
   virtual fastNLOCoeffBase* Clone() const;                                     //!< returns 'new' copy of this instance.
   virtual void Read(istream& table);
   virtual void Write(ostream& table);
   virtual void Print() const;
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false);

protected:
   fastNLOCoeffData();
   void ReadCoeffData(istream& table);
   void ReadRest(istream& table);

   int Nuncorrel;
   vector<string > UncDescr;
   int Ncorrel;
   vector<string > CorDescr;
   vector<double > Xcenter;
   vector<double > Value;
   v2d UncorLo;
   v2d UncorHi;
   v2d CorrLo;
   v2d CorrHi;
   int NErrMatrix;
   v2d matrixelement;

};

#endif
