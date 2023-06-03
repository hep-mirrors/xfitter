#ifndef __fastNLOCoeffMult__
#define __fastNLOCoeffMult__

#include "fastNLOCoeffBase.h"
#include "fastNLOConstants.h"


class fastNLOCoeffMult : public fastNLOCoeffBase {

   friend class fastNLOTable;
   friend class fastNLOCreate;

public:
   fastNLOCoeffMult() = delete;
   fastNLOCoeffMult(int NObsBin);
   explicit fastNLOCoeffMult(const fastNLOCoeffBase&);
   virtual ~fastNLOCoeffMult(){;};
   virtual fastNLOCoeffMult* Clone() const;                                     //!< returns 'new' copy of this instance.
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false);
   virtual void Read(std::istream& table, int ITabVersionRead);
   virtual void Write(std::ostream& table, int ITabVersionWrite);
   virtual void Print(int iprint) const;

   double GetMultFactor(int iObs) const { return fact[iObs]; }
   std::vector<double > GetMultFactor() const { return fact; }
   std::vector<std::string> GetUncDescription() const { return UncDescr; }
   std::vector<std::string> GetCorDescription() const { return CorDescr; }
   fastNLO::v2d GetUncorLo() const { return UncorLo; };
   fastNLO::v2d GetUncorHi() const { return UncorHi; };
   fastNLO::v2d GetCorrLo()  const { return CorrLo; };
   fastNLO::v2d GetCorrHi()  const { return CorrHi; };

   // Erase observable bin; iObsIdx is the C++ array index to be removed and
   // not the observable bin no. running from 1 to fNObsBins
   virtual void EraseBin(unsigned int iObsIdx);
   virtual void MultiplyBin(unsigned int iObsIdx, double fact);
   bool IsCatenable(const fastNLOCoeffMult& other) const;
   // Catenate observable to table
   virtual void CatBin(const fastNLOCoeffMult& other, unsigned int iObsIdx);

   // getter/setter
   int  GetNuncorrel() const {return Nuncorrel;}
   void SetNuncorrel(int n){Nuncorrel = n;}
   int  GetNcorrel() const {return Ncorrel;}
   void SetNcorrel(int n){Ncorrel = n;}

protected:
   void ReadCoeffMult(std::istream& table);
   void ReadRest(std::istream& table, int ITabVersionRead);

   int Nuncorrel;
   std::vector < std::string > UncDescr;
   int Ncorrel;
   std::vector < std::string > CorDescr;
   fastNLO::v2d UncorLo;
   fastNLO::v2d UncorHi;
   fastNLO::v2d CorrLo;
   fastNLO::v2d CorrHi;
   fastNLO::v1d fact;

};

#endif
