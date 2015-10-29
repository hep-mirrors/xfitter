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

using namespace std;
using namespace fastNLO;

class fastNLOCoeffBase : public PrimalScream {

   friend class fastNLOTable;

public:
   fastNLOCoeffBase(int NObsBin);                                               //! Use this constructor
   virtual ~fastNLOCoeffBase(){;};                                              //! destructor
   //fastNLOCoeffBase(const fastNLOCoeffBase& coeff);                           //! Use compiler-default
   virtual fastNLOCoeffBase* Clone() const;                                     //!< returns 'new' copy of this instance.

   virtual void Read(istream& table);
   virtual void Write(ostream& table);
   //void Add(fastNLOCoeffBase* other);
   virtual void Print() const;

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

   bool GetIsFlexibleScale() const { return (NScaleDep>=3) && (IAddMultFlag==0); }

   vector<string > GetContributionDescription() const { return CtrbDescript; }
   void SetContributionDescription(vector<string > descr ) { CtrbDescript = descr; };           //! Set contribution description
   vector<string > GetCodeDescription() const { return CodeDescript; }


   bool IsLO() const {return IContrFlag1==1 && IContrFlag2==1;}
   bool IsNLO() const {return IContrFlag1==1 && IContrFlag2==2;}
   bool IsNNLO() const {return IContrFlag1==1 && IContrFlag2==3;}

   bool IsCompatible(const fastNLOCoeffBase& other) const;
   //bool operator==(const fastNLOCoeffBase& other) const { return IsCompatible(other); }


protected:
   fastNLOCoeffBase();
   void ReadBase(istream& table);
   void EndReadCoeff(istream& table);

   int fNObsBins; // obtained from Scenario

   int IXsectUnits;
   int IDataFlag;
   int IAddMultFlag;
   int IContrFlag1;
   int IContrFlag2;
   int NScaleDep;
   vector < string > CtrbDescript;
   vector < string > CodeDescript;

};


#endif
