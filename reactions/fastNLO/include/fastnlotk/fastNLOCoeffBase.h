#ifndef __fastNLOCoeffBase__
#define __fastNLOCoeffBase__

#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <iostream>

#include "fastNLOConstants.h"
#include "speaker.h"


class fastNLOCoeffBase : public PrimalScream {

   friend class fastNLOTable;

public:
   fastNLOCoeffBase() = delete;
   fastNLOCoeffBase(int NObsBin);                                               //! Use this constructor
   // deletes instance of derived classes through pointer to base class
   virtual ~fastNLOCoeffBase(){};                                              //! destructor
   virtual fastNLOCoeffBase* Clone() const;                                     //!< returns 'new' copy of this instance.

   virtual void Read(std::istream& table);
   virtual void Write(std::ostream& table, int ITabVersionWrite);
   virtual void Print(int iprint) const;

   // Erase or multiply observable bin; iObsIdx is the C++ array index to be removed and
   // not the observable bin no. running from 1 to fNObsBins
   virtual void EraseBin(unsigned int iObsIdx);
   virtual void MultiplyBin(unsigned int iObsIdx, double fact);
   // Catenate observable to table
   virtual void CatBin(const fastNLOCoeffBase& other, unsigned int iObsIdx);

   bool IsCatenable(const fastNLOCoeffBase& other) const;

   void SetCoeffAddDefaults();

   int GetIDataFlag() const {return IDataFlag;}
   void SetIDataFlag(int n){IDataFlag = n;}

   int GetIAddMultFlag() const {return IAddMultFlag;}
   void SetIAddMultFlag(int n){IAddMultFlag = n;}

   int GetIContrFlag1() const {return IContrFlag1;}
   void SetIContrFlag1(int n){IContrFlag1 = n;}

   int GetIContrFlag2() const {return IContrFlag2;}
   void SetIContrFlag2(int n){IContrFlag2 = n;}

   int GetNScaleDep() const {return NScaleDep;}
   void SetNScaleDep(int n){NScaleDep = n;}

   int GetIXsectUnits() const { return IXsectUnits;}
   void SetIXsectUnits(int n){IXsectUnits = n;}

   int GetNObsBin() const { return fNObsBins;}
   void SetNObsBin(unsigned int nObs) { fNObsBins = nObs;}

   bool GetIsFlexibleScale() const { return (NScaleDep>=3) && (IAddMultFlag==0); }

   std::vector<std::string > GetContributionDescription() const { return CtrbDescript; }
   void SetContributionDescription(std::vector<std::string > descr ) { CtrbDescript = descr; };           //! Set contribution description
   std::vector<std::string > GetCodeDescription() const { return CodeDescript; }

   bool IsLO() const {return IContrFlag1==1 && IContrFlag2==1;}
   bool IsNLO() const {return IContrFlag1==1 && IContrFlag2==2;}
   bool IsNNLO() const {return IContrFlag1==1 && IContrFlag2==3;}
   bool IsCompatible(const fastNLOCoeffBase& other) const;

   bool IsEnabled() const {return enabled;}
   void Enable(bool on=true) {enabled = on;}

protected:
   void ReadBase(std::istream& table);
   void EndReadCoeff(std::istream& table);

   int fNObsBins; // obtained from Scenario

   int IXsectUnits;
   int IDataFlag;
   int IAddMultFlag;
   int IContrFlag1;
   int IContrFlag2;
   int NScaleDep;
   int fVersionRead = 23000;
   std::vector < std::string > CtrbDescript;
   std::vector < std::string > CodeDescript;

   bool enabled = false;
};


#endif
