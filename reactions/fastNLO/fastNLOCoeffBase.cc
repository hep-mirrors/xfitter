#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOCoeffBase.h"
#include "fastnlotk/fastNLOConstants.h"
#include "fastnlotk/fastNLOTools.h"

using namespace std;
using namespace fastNLO;

//________________________________________________________________________________________________________________ //
fastNLOCoeffBase::fastNLOCoeffBase(int NObsBin)
   : PrimalScream("fastNLOCoeffBase"), fNObsBins(NObsBin), IXsectUnits(),
     IDataFlag(), IAddMultFlag(), IContrFlag1(), IContrFlag2(), NScaleDep(),
     CtrbDescript(), CodeDescript() {}

//________________________________________________________________________________________________________________ //
fastNLOCoeffBase* fastNLOCoeffBase::Clone() const {
   //! User has to take care to delete this object later
   return new fastNLOCoeffBase(*this);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::Read(istream& table, int ITabVersionRead){
   // basic read function.
   // reads in only 'base'-variables.
   debug["Read"]<<"Start reading table ..."<<endl;
   ReadBase(table, ITabVersionRead);
   ReadCoeffInfoBlocks(table, ITabVersionRead);
   EndReadCoeff(table, ITabVersionRead);
   debug["Read"]<<"Finished reading table."<<endl;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::ReadBase(istream& table, int ITabVersionRead){
   debug["ReadBase"]<<"Start reading base coefficient table ..."<<endl;
   fastNLOTools::ReadMagicNo(table);
   table >> IXsectUnits;
   table >> IDataFlag;
   table >> IAddMultFlag;
   table >> IContrFlag1;
   table >> IContrFlag2;
   table >> NScaleDep;
   fastNLOTools::ReadFlexibleVector(CtrbDescript,table);
   fastNLOTools::ReadFlexibleVector(CodeDescript,table);
   debug["ReadBase"]<<"Finished reading base coefficient table."<<endl;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::EndReadCoeff(istream& table, int ITabVersionRead){
   debug["EndReadCoeff"]<<"Should have reached end of coefficient table for table version "<<ITabVersionRead<<endl;
   fastNLOTools::ReadMagicNo(table);
   fastNLOTools::PutBackMagicNo(table);
   debug["EndReadCoeff"]<<"Finished reading coefficient table for table version "<<ITabVersionRead<<endl;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::Write(ostream& table, int ITabVersionWrite) {
   say::debug["Write"]<<"Writing fastNLOCoeffBase for table version " << ITabVersionWrite << "." << endl;
   table << fastNLO::tablemagicno << sep;
   // if ( itabversion >= 24000 ) table << fastNLO::tabversion << sep;
   // if ( itabversion >= 24000 ) table << "fastNLO_CoeffAddBase" << sep;
   // if ( itabversion >= 24000 ) table << 0 << sep; // v2.4, but yet unused
   table << IXsectUnits << sep;
   table << IDataFlag << sep;
   table << IAddMultFlag << sep;
   table << IContrFlag1 << sep;
   table << IContrFlag2 << sep;
   table << NScaleDep << sep;
   fastNLOTools::WriteFlexibleVector(CtrbDescript,table);
   fastNLOTools::WriteFlexibleVector(CodeDescript,table);
   // if ( itabversion >= 24000 ) table << 0 << sep; // v2.4, but yet unuse
}


//________________________________________________________________________________________________________________ //
bool fastNLOCoeffBase::IsCompatible(const fastNLOCoeffBase& other) const {
   if( fNObsBins != other.GetNObsBin() ){
      warn["IsCompatible"]<<"fNObsBins != other.GetNObsBin()"<<endl;
      return false;
   }
   if( IXsectUnits != other.GetIXsectUnits() ){
      warn["IsCompatible"]<<"IXsectUnits != other.GetIXsectUnits()"<<endl;
      return false;
   }
   if( IDataFlag != other.GetIDataFlag() ){
      debug["IsCompatible"]<<"IDataFlag != other.GetIDataFlag()"<<endl;
      return false;
   }
   if( IAddMultFlag != other.GetIAddMultFlag() ){
      debug["IsCompatible"]<<"IAddMultFlag != other.GetIAddMultFlag()"<<endl;
      return false;
   }
   if( IContrFlag1 != other.GetIContrFlag1() ){
      debug["IsCompatible"]<<"IContrFlag1 != other.GetIContrFlag1()"<<endl;
      return false;
   }
   if( IContrFlag2 != other.GetIContrFlag2() ){
      debug["IsCompatible"]<<"IContrFlag2 != other.GetIContrFlag2()"<<endl;
      return false;
   }
   if( NScaleDep != other.GetNScaleDep() ){
      debug["IsCompatible"]<<"NScaleDep != other.GetNScaleDep()"<<endl;
      if ( (NScaleDep==5 && other.GetNScaleDep()==6) || (NScaleDep==6 && other.GetNScaleDep()==5) ) {
         debug["IsCompatible"]<<"One table with NScale=5 and one with NScaleDep=6"<<endl;
         // continue;
      }
      else {
         warn["IsCompatible"]<<"Incompatible NScaleDep found!()"<<endl;
         return false;
      }
   }
   debug["IsCompatible"]<<"Both tables are compatible"<<endl;
   // check descripts here ?!
   //bool potentialcompatible = true;
   //vector < string > CtrbDescript;
   //vector < string > CodeDescript;
   return true;
}


//________________________________________________________________________________________________________________ //
bool fastNLOCoeffBase::IsCatenable(const fastNLOCoeffBase& other) const {
   if ( IXsectUnits != other.GetIXsectUnits() ){
      debug["IsCatenable"]<<"IXsectUnits != other.GetIXsectUnits(). Skipped."<<endl;
      return false;
   }
   if ( IDataFlag != other.GetIDataFlag() ){
      debug["IsCatenable"]<<"IDataFlag != other.GetIDataFlag(). Skipped."<<endl;
      return false;
   }
   if ( IAddMultFlag != other.GetIAddMultFlag() ){
      debug["IsCatenable"]<<"IAddMultFlag != other.GetIAddMultFlag(). Skipped."<<endl;
      return false;
   }
   if ( IContrFlag1 != other.GetIContrFlag1() ){
      debug["IsCatenable"]<<"IContrFlag1 != other.GetIContrFlag1(). Skipped."<<endl;
      return false;
   }
   if ( IContrFlag2 != other.GetIContrFlag2() ){
      debug["IsCatenable"]<<"IContrFlag2 != other.GetIContrFlag2(). Skipped."<<endl;
      return false;
   }
   if ( NScaleDep != other.GetNScaleDep() ){
      return false;
   }
   if ( ! (fVersionRead < 25000) ) {
      if ( (HasCoeffInfoBlock(0,0) && ! other.HasCoeffInfoBlock(0,0)) ||
           (! HasCoeffInfoBlock(0,0) && other.HasCoeffInfoBlock(0,0)) ) {
         debug["IsCatenable"]<<"Missing InfoBlock in either of the two tables. Skipped."<<endl;
         return false;
      }
   }
   info["IsCatenable"]<<"Base parameters of contribution allow catenation"<<endl;
   // check descripts here ?!
   //bool potentialcompatible = true;
   //vector < string > CtrbDescript;
   //vector < string > CodeDescript;
   return true;
}

//________________________________________________________________________________________________________________ //


void fastNLOCoeffBase::SetCoeffAddDefaults(){
  SetIDataFlag(0);
  SetIAddMultFlag(0);
  SetIContrFlag1(1);
  SetIContrFlag2(100); // specify if LO or NLO
  SetNScaleDep(0);
  SetIXsectUnits(12);
};

//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::Print(int iprint) const {
   if ( !(iprint < 0) ) {
      cout << fastNLO::_DSEP20C << " fastNLO Table: CoeffBase " << fastNLO::_DSEP20 << endl;
   } else {
      cout << endl << fastNLO::_CSEP20C << " fastNLO Table: CoeffBase " << fastNLO::_CSEP20 << endl;
   }
   fastNLOTools::PrintVector(CtrbDescript,"Contribution description (CtrbDescript)","#");
   fastNLOTools::PrintVector(CodeDescript,"Code description (CodeDescript)","#");
   if ( std::abs(iprint) > 0 ) {
      cout << fastNLO::_SSEP20C << " Extended information (iprint > 0) " << fastNLO::_SSEP20 << endl;
      printf(" #   IXsectUnits                       %d\n",IXsectUnits);
      printf(" #   IDataFlag                         %d\n",IDataFlag);
      printf(" #   IAddMultFlag                      %d\n",IAddMultFlag);
      printf(" #   IContrFlag1                       %d\n",IContrFlag1);
      printf(" #   IContrFlag2                       %d\n",IContrFlag2);
      printf(" #   NScaleDep                         %d\n",NScaleDep);
   }
   if ( iprint < 0 ) {
      cout << fastNLO::_CSEPSC << endl;
   } else {
      //      cout << fastNLO::_DSEPSC << endl;
   }
}


//________________________________________________________________________________________________________________ //

// Erase observable bin
void fastNLOCoeffBase::EraseBin(unsigned int iObsIdx) {
   debug["EraseBin"]<<"Erasing table entries in CoeffBase for bin index " << iObsIdx << endl;
   SetNObsBin(GetNObsBin()-1);
}

// Catenate observable bin
void fastNLOCoeffBase::CatBin(const fastNLOCoeffBase& other, unsigned int iObsIdx) {
   debug["CatBin"]<<"Catenating observable bin in CoeffBase corresponding to bin index " << iObsIdx << endl;
   SetNObsBin(GetNObsBin()+1);
}

// Multiply observable bin
void fastNLOCoeffBase::MultiplyBin(unsigned int iObsIdx, double nfact) {
   debug["MultiplyBin"]<<"Multiplying table entries. Nothing to be done in CoeffBase." << endl;
}

//
//________________________________________________________________________________________________________________
// Added to include CoeffInfoBlocks
// Check existence of InfoBlock of type, i.e. flags, (flag1,x); return error in case of multiple matches!
bool fastNLOCoeffBase::HasCoeffInfoBlock(int fICoeffInfoBlockFlag1) const {
   bool result = false;
   for (int i=0; i<NCoeffInfoBlocks; i++) {
      if ( ICoeffInfoBlockFlag1[i] == fICoeffInfoBlockFlag1 ) {
         if ( ! result ) {
            result = true;
         } else {
            error["HasCoeffInfoBlocks"]<<"Aborted! Found multiple info blocks of type ICoeffInfoBlockFlag1 = "<<ICoeffInfoBlockFlag1[i]<<endl;
            exit(201);
         }
      }
   }
   return result;
}

// Check existence of InfoBlock of type, i.e. flags, (flag1,flag2); return error in case of multiple matches!
bool fastNLOCoeffBase::HasCoeffInfoBlock(int fICoeffInfoBlockFlag1, int fICoeffInfoBlockFlag2) const {
   bool result = false;
   for (int i=0; i<NCoeffInfoBlocks; i++) {
      if ( ICoeffInfoBlockFlag1[i] == fICoeffInfoBlockFlag1 && ICoeffInfoBlockFlag2[i] == fICoeffInfoBlockFlag2 ) {
         if ( ! result ) {
            result = true;
         } else {
            error["HasCoeffInfoBlocks"]<<"Aborted! Found multiple info blocks of type ICoeffInfoBlockFlag1 = "<<ICoeffInfoBlockFlag1[i]<<endl;
            exit(202);
         }
      }
   }
   return result;
}

// Return first matching index; multiple blocks of same type, i.e. flags (flag1,x), are not allowed!
int fastNLOCoeffBase::GetCoeffInfoBlockIndex(int fICoeffInfoBlockFlag1) {
   for (int i=0; i<NCoeffInfoBlocks; i++) {
      if ( ICoeffInfoBlockFlag1[i] == fICoeffInfoBlockFlag1 ) return i;
   }
   return -1;
}

// Return first matching index; multiple blocks of same type, i.e. flags (flag1,flag2), are not allowed!
int fastNLOCoeffBase::GetCoeffInfoBlockIndex(int fICoeffInfoBlockFlag1, int fICoeffInfoBlockFlag2) {
   for (int i=0; i<NCoeffInfoBlocks; i++) {
      if ( ICoeffInfoBlockFlag1[i] == fICoeffInfoBlockFlag1 && ICoeffInfoBlockFlag2[i] == fICoeffInfoBlockFlag2 ) return i;
   }
   return -1;
}

void fastNLOCoeffBase::AddCoeffInfoBlock(int fICoeffInfoBlockFlag1, int fICoeffInfoBlockFlag2, std::vector<std::string> Description,
                                         std::vector<double> Uncertainty) {
   info["AddCoeffInfoBlocks"]<<"Adding additional InfoBlock with flags "<<fICoeffInfoBlockFlag1<<" and "<<fICoeffInfoBlockFlag2<<" to table contribution."<<endl;
   NCoeffInfoBlocks += 1;
   ICoeffInfoBlockFlag1.push_back(fICoeffInfoBlockFlag1);
   ICoeffInfoBlockFlag2.push_back(fICoeffInfoBlockFlag2);
   NCoeffInfoBlockDescr.push_back(Description.size());
   CoeffInfoBlockDescript.push_back(Description);
   NCoeffInfoBlockCont.push_back(Uncertainty.size());
   CoeffInfoBlockContent.push_back(Uncertainty);
}

void fastNLOCoeffBase::AddCoeffInfoBlock(int fICoeffInfoBlockFlag1, int fICoeffInfoBlockFlag2, std::vector<std::string> Description,
                                         std::string filename, unsigned int icola, unsigned int icolb) {
   info["AddCoeffInfoBlocks"]<<"Adding additional InfoBlock reading data from file "<<filename<<endl;
   std::vector<double> Uncertainty = fastNLOTools::ReadUncertaintyFromFile(filename, icola, icolb);
   AddCoeffInfoBlock(fICoeffInfoBlockFlag1, fICoeffInfoBlockFlag2, Description, Uncertainty);
}

void fastNLOCoeffBase::ReadCoeffInfoBlocks(istream& table, int ITabVersionRead) {
   if (ITabVersionRead < 25000) {
      debug["ReadCoeffInfoBlocks"]<<"No additional info blocks allowed for table versions < 25000"<<endl;
   } else {
      table >> NCoeffInfoBlocks;
      debug["ReadCoeffInfoBlocks"]<<"Found "<<NCoeffInfoBlocks<<" additional info blocks for coefficient table."<<endl;
      for (int i=0; i<NCoeffInfoBlocks; i++) {
         int iflag;
         table >> iflag;
         ICoeffInfoBlockFlag1.push_back(iflag);
         if (ICoeffInfoBlockFlag1[i] == 0) {
            debug["ReadCoeffInfoBlocks"]<<"Found info block of type ICoeffInfoBlockFlag1 = "<<ICoeffInfoBlockFlag1[i]<<endl;
         } else {
            error["ReadCoeffInfoBlocks"]<<"Found info block of unknown type ICoeffInfoBlockFlag1 = "<<ICoeffInfoBlockFlag1[i]<<endl;
            exit(111);
         }
         table >> iflag;
         ICoeffInfoBlockFlag2.push_back(iflag);
         if (ICoeffInfoBlockFlag2[i] == 0 || ICoeffInfoBlockFlag2[i] == 1) {
            debug["ReadCoeffInfoBlocks"]<<"Found info block of type ICoeffInfoBlockFlag2 = "<<ICoeffInfoBlockFlag2[i]<<endl;
         } else {
            error["ReadCoeffInfoBlocks"]<<"Found info block of unknown type ICoeffInfoBlockFlag2 = "<<ICoeffInfoBlockFlag2[i]<<endl;
            exit(222);
         }
         std::vector < std::string > Description;
         NCoeffInfoBlockDescr.push_back(fastNLOTools::ReadFlexibleVector(Description,table)-1);
         CoeffInfoBlockDescript.push_back(Description);
         for (unsigned int j=0; j<Description.size(); j++) {
            debug["ReadCoeffInfoBlocks"]<<"Read info block description line "<<j<<" : "<< Description[j] <<endl;
         }
         if ( ICoeffInfoBlockFlag1[i] == 0 ) { // Entry per ObsBin
            std::vector < double > Content;
            NCoeffInfoBlockCont.push_back(fastNLOTools::ReadFlexibleVector(Content,table));
            if ( NCoeffInfoBlockCont[i]-1 != fNObsBins ) {
               error["ReadCoeffInfoBlocks"]<<"Found info block of type ICoeffInfoBlockFlag1 = "<<ICoeffInfoBlockFlag1[i]<<" , but # of content lines = "<<NCoeffInfoBlockCont[i]-1<<" differs from fNObsBins = "<<fNObsBins<<"! Aborted."<<endl;
               exit(223);
            } else {
               debug["ReadCoeffInfoBlocks"]<<"Read "<<NCoeffInfoBlockCont[i]-1<<" lines into InfoBlock content vector."<<endl;
            }
            CoeffInfoBlockContent.push_back(Content);
         }
      }
   }
}

void fastNLOCoeffBase::WriteCoeffInfoBlocks(ostream& table, int ITabVersionWrite) {
   if (ITabVersionWrite < 25000) {
      debug["WriteCoeffInfoBlocks"]<<"No additional InfoBlocks allowed for table versions < 25000"<<endl;
   } else {
      debug["WriteCoeffInfoBlocks"]<<"Writing additional InfoBlocks; NCoeffInfoBlocks = "<<NCoeffInfoBlocks<<endl;
      table << NCoeffInfoBlocks << sep;
      for (int i=0; i<NCoeffInfoBlocks; i++) {
         debug["WriteCoeffInfoBlocks"]<<"ICoeffInfoBlockFlags1,2 = "<<ICoeffInfoBlockFlag1[i]<<", "<<ICoeffInfoBlockFlag1[i]<<endl;
         table << ICoeffInfoBlockFlag1[i] << sep;
         table << ICoeffInfoBlockFlag2[i] << sep;
         table << NCoeffInfoBlockDescr[i] << sep;
         for (int j=0; j<NCoeffInfoBlockDescr[i]; j++) {
            debug["WriteCoeffInfoBlocks"]<<"CoeffInfoBlockDecsript[.][.] = "<<CoeffInfoBlockDescript[i][j]<<endl;
            table << CoeffInfoBlockDescript[i][j] << sep;
         }
         int nlines = fastNLOTools::WriteFlexibleVector(CoeffInfoBlockContent[i],table);
         debug["WriteCoeffInfoBlocks"]<<"Wrote "<<nlines-1<<" content lines to additional InfoBlock no. "<<i<<endl;
      }
   }
}

void fastNLOCoeffBase::SetContributionDescription(std::vector <std::string> CtrbDescr) {
   debug["SetContributionDescription"]<<"Setting contribution description."<<endl;
   size_t NCtrbDescript = CtrbDescr.size();
   fastNLOCoeffBase::CtrbDescript.resize(NCtrbDescript);
   for (size_t i=0; i < NCtrbDescript; ++i) {
      fastNLOCoeffBase::CtrbDescript[i] = CtrbDescr[i];
   }
}

void fastNLOCoeffBase::SetCodeDescription(std::vector <std::string> CodeDescr) {
   debug["SetCodeDescription"]<<"Setting code description."<<endl;
   size_t NCodeDescript = CodeDescr.size();
   fastNLOCoeffBase::CodeDescript.resize(NCodeDescript);
   for (size_t i=0; i < NCodeDescript; ++i) {
      fastNLOCoeffBase::CodeDescript[i] = CodeDescr[i];
   }
}
