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

   virtual void Read(std::istream& table, int ITabVersionRead);
   virtual void Write(std::ostream& table, int ITabVersionWrite);
   virtual void Print(int iprint) const;

   // Erase or multiply observable bin; iObsIdx is the C++ array index to be removed and
   // not the observable bin no. running from 1 to fNObsBins
   virtual void EraseBin(unsigned int iObsIdx);
   virtual void MultiplyBin(unsigned int iObsIdx, double fact);
   // Catenate observable to table
   virtual void CatBin(const fastNLOCoeffBase& other, unsigned int iObsIdx);

   bool IsCatenable(const fastNLOCoeffBase& other) const;
   virtual bool IsEquivalent(const fastNLOCoeffBase& other, double rtol) const;

   void SetCoeffAddDefaults();

   /// ___________________________________________________________________________________________________
   /// Some info getters & setters for contribution modifications
   /// ___________________________________________________________________________________________________

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

   /// get/set contribution and code description
   std::vector<std::string > GetContributionDescription() const { return CtrbDescript; }
   void SetContributionDescription(std::vector<std::string > descr );
   std::vector<std::string > GetCodeDescription() const { return CodeDescript; }
   void SetCodeDescription(std::vector<std::string > descr );

   bool IsLO() const {return IContrFlag1==1 && IContrFlag2==1;}
   bool IsNLO() const {return IContrFlag1==1 && IContrFlag2==2;}
   bool IsNNLO() const {return IContrFlag1==1 && IContrFlag2==3;}
   bool IsCompatible(const fastNLOCoeffBase& other) const;

   bool IsEnabled() const {return enabled;}
   void Enable(bool on=true) {enabled = on;}

   // Added to include CoeffInfoBlocks
   bool HasCoeffInfoBlock() const {return NCoeffInfoBlocks>0;}
   bool HasCoeffInfoBlock(int ICoeffInfoBlockFlag1) const;
   bool HasCoeffInfoBlock(int ICoeffInfoBlockFlag1, int ICoeffInfoBlockFlag2) const;
   int GetCoeffInfoBlockIndex(int ICoeffInfoBlockFlag1);
   int GetCoeffInfoBlockIndex(int ICoeffInfoBlockFlag1, int ICoeffInfoBlockFlag2);
   int GetCoeffInfoBlockFlag1(int Index) const { return ICoeffInfoBlockFlag1[Index]; };
   int GetCoeffInfoBlockFlag2(int Index) const { return ICoeffInfoBlockFlag2[Index]; };
   void SetCoeffInfoBlockFlag1(int Index, int iFlag1) { ICoeffInfoBlockFlag1[Index] = iFlag1; };
   void SetCoeffInfoBlockFlag2(int Index, int iFlag2) { ICoeffInfoBlockFlag2[Index] = iFlag2; };
   std::vector< std::string > GetCoeffInfoBlockDescription(int Index) const { return CoeffInfoBlockDescript[Index]; }
   std::vector < double > GetCoeffInfoBlockContent(int Index) const { return CoeffInfoBlockContent[Index]; };
   int GetNCoeffInfoBlocks() const {return NCoeffInfoBlocks;}
   // Provide uncertainty via input vector
   void AddCoeffInfoBlock(int ICoeffInfoBlockFlag1, int ICoeffInfoBlockFlag2, std::vector<std::string> Description,
                          std::vector<double> Content);
   // Provide uncertainty reading from filename
   void AddCoeffInfoBlock(int ICoeffInfoBlockFlag1, int ICoeffInfoBlockFlag2, std::vector<std::string> Description,
                          std::string filename, unsigned int icola = 0, unsigned int icolb = 0, double relfac = 1);



protected:
   void ReadBase(std::istream& table, int ITabVersionRead);
   void ReadCoeffInfoBlocks(std::istream& table, int ITabVersionRead);
   void WriteCoeffInfoBlocks(std::ostream& table, int ITabVersionWrite);
   void EndReadCoeff(std::istream& table, int ITabVersionRead);

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

   // Added to include CoeffInfoBlocks
   int NCoeffInfoBlocks = 0; // Not present for version numbers < 25000
   std::vector < int > ICoeffInfoBlockFlag1;
   std::vector < int > ICoeffInfoBlockFlag2;
   std::vector < int > NCoeffInfoBlockDescr;
   std::vector < std::vector < std::string > > CoeffInfoBlockDescript;
   std::vector < int > NCoeffInfoBlockCont;
   fastNLO::v2d CoeffInfoBlockContent;
};


#endif
