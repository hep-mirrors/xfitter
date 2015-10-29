#ifndef __fastNLOCoeffMult__
#define __fastNLOCoeffMult__

#include "fastNLOCoeffBase.h"
#include "fastNLOConstants.h"

using namespace std;

class fastNLOCoeffMult : public fastNLOCoeffBase {

   friend class fastNLOTable;
   friend class fastNLOCreate;

public:
   fastNLOCoeffMult();
   fastNLOCoeffMult(int NObsBin);
   fastNLOCoeffMult(const fastNLOCoeffBase&);
   virtual ~fastNLOCoeffMult(){;};
   virtual fastNLOCoeffBase* Clone() const;                                     //!< returns 'new' copy of this instance.
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false);
   virtual void Read(istream& table);
   virtual void Write(ostream& table);
   virtual void Print() const;

   double GetMultFactor(int iObs) const { return fact[iObs]; }
   vector<double > GetMultFactor() const { return fact; }
   vector<string> GetUncDescription() const { return UncDescr; }
   vector<string> GetCorDescription() const { return CorDescr; }
   v2d GetUncorLo() const { return UncorHi; };
   v2d GetUncorHi() const { return UncorLo; };
   v2d GetCorrLo()  const { return CorrLo; };
   v2d GetCorrHi()  const { return CorrHi; };

protected:
   void ReadCoeffMult(istream& table);
   void ReadRest(istream& table);

   int Nuncorrel;
   vector < string > UncDescr;
   int Ncorrel;
   vector < string > CorDescr;
   v2d UncorLo;
   v2d UncorHi;
   v2d CorrLo;
   v2d CorrHi;
   v1d fact;

};

#endif
