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
void fastNLOCoeffBase::Read(istream& table){
   // basic read function.
   // reads in only 'base'-variables.
   debug["Read"]<<endl;
   ReadBase(table);
   EndReadCoeff(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::ReadBase(istream& table){
   debug["ReadBase"]<<endl;
   //table.peek();

   table >> fVersionRead;
   if ( fVersionRead == fastNLO::tablemagicno )
      fVersionRead = 22000 ;
   //fastNLOTools::ReadMagicNo(table);
   std::string stest;
   if ( fVersionRead>=24000 ) table >> stest; //"fastNLO_CoeffAddBase"
   if ( fVersionRead>=24000 ) fastNLOTools::ReadUnused(table);

   table >> IXsectUnits;
   table >> IDataFlag;
   table >> IAddMultFlag;
   table >> IContrFlag1;
   table >> IContrFlag2;
   table >> NScaleDep;
   fastNLOTools::ReadFlexibleVector(CtrbDescript,table);
   fastNLOTools::ReadFlexibleVector(CodeDescript,table);
   //printf("  *  fastNLOCoeffBase::Read().  IDataFlag: %d, IAddMultFlag: %d, IContrFlag1: %d, IContrFlag2: %d,, NScaleDep: %d\n",IDataFlag,IAddMultFlag,IContrFlag1,IContrFlag2,NScaleDep );
   // int NContrDescr;
   // table >> NContrDescr;
   // CtrbDescript.resize(NContrDescr);
   // char buffer[5257];
   // table.getline(buffer,5256);
   // for(int i=0;i<NContrDescr;i++){
   //    table.getline(buffer,256);
   //    CtrbDescript[i] = buffer;
   //    //      StripWhitespace(CtrbDescript[i]);
   // }
   // int NCodeDescr;
   // table >> NCodeDescr;
   // CodeDescript.resize(NCodeDescr);
   // table.getline(buffer,256);
   // for(int i=0;i<NCodeDescr;i++){
   //    table.getline(buffer,256);
   //    CodeDescript[i] = buffer;
   //    //      StripWhitespace(CodeDescript[i]);
   // }
   if ( fVersionRead>=24000 ) fastNLOTools::ReadUnused(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::EndReadCoeff(istream& table){
   fastNLOTools::ReadMagicNo(table);
   // if (!fastNLOTools::ReadMagicNo(table)) {
   //    say::error["ReadBase"]<<"Did not find final magic number, aborting!"<<endl;
   //    say::error["ReadBase"]<<"Please check compatibility of tables and program version!"<<endl;
   //    say::error["ReadBase"]<<"This might also be provoked by lines with unexpected non-numeric content like 'inf' or 'nan'!"<<endl;
   //    exit(1);
   // }
   fastNLOTools::PutBackMagicNo(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::Write(ostream& table, int itabversion) {
   say::debug["Write"]<<"Writing fastNLOCoeffBase for table version " << itabversion << "." << endl;
   table << fastNLO::tablemagicno << sep;
   if ( itabversion >= 24000 ) table << fastNLO::tabversion << sep;
   if ( itabversion >= 24000 ) table << "fastNLO_CoeffAddBase" << sep;
   if ( itabversion >= 24000 ) table << 0 << sep; // v2.4, but yet unused
   table << IXsectUnits << sep;
   table << IDataFlag << sep;
   table << IAddMultFlag << sep;
   table << IContrFlag1 << sep;
   table << IContrFlag2 << sep;
   table << NScaleDep << sep;
   //printf("  *  fastNLOCoeffBase::Write().  IDataFlag: %d, IAddMultFlag: %d, IContrFlag1: %d, IContrFlag2: %d, NScaleDep: %d\n",
   //IDataFlag,IAddMultFlag,IContrFlag1,IContrFlag2,NScaleDep);
   fastNLOTools::WriteFlexibleVector(CtrbDescript,table);
   fastNLOTools::WriteFlexibleVector(CodeDescript,table);
   if ( itabversion >= 24000 ) table << 0 << sep; // v2.4, but yet unuse
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
   if( IXsectUnits != other.GetIXsectUnits() ){
      debug["IsCatenable"]<<"IXsectUnits != other.GetIXsectUnits(). Skipped."<<endl;
      return false;
   }
   if( IDataFlag != other.GetIDataFlag() ){
      debug["IsCatenable"]<<"IDataFlag != other.GetIDataFlag(). Skipped."<<endl;
      return false;
   }
   if( IAddMultFlag != other.GetIAddMultFlag() ){
      debug["IsCatenable"]<<"IAddMultFlag != other.GetIAddMultFlag(). Skipped."<<endl;
      return false;
   }
   if( IContrFlag1 != other.GetIContrFlag1() ){
      debug["IsCatenable"]<<"IContrFlag1 != other.GetIContrFlag1(). Skipped."<<endl;
      return false;
   }
   if( IContrFlag2 != other.GetIContrFlag2() ){
      debug["IsCatenable"]<<"IContrFlag2 != other.GetIContrFlag2(). Skipped."<<endl;
      return false;
   }
   if( NScaleDep != other.GetNScaleDep() ){
      debug["IsCatenable"]<<"NScaleDep != other.GetNScaleDep(). Skipped."<<endl;
      return false;
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
   if ( abs(iprint) > 0 ) {
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
   debug["fastNLOCoeffBase::EraseBin"]<<"Erasing table entries in CoeffBase for bin index " << iObsIdx << endl;
   SetNObsBin(GetNObsBin()-1);
}

// Catenate observable bin
void fastNLOCoeffBase::CatBin(const fastNLOCoeffBase& other, unsigned int iObsIdx) {
   debug["fastNLOCoeffBase::CatBin"]<<"Catenating observable bin in CoeffBase corresponding to bin index " << iObsIdx << endl;
   SetNObsBin(GetNObsBin()+1);
}

// Multiply observable bin
void fastNLOCoeffBase::MultiplyBin(unsigned int iObsIdx, double nfact) {
   debug["fastNLOCoeffBase::MultiplyBin"]<<"Multiplying table entries. Nothing to be done in CoeffBase." << endl;
}
