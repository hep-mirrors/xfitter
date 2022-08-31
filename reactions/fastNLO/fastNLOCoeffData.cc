#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOCoeffData.h"
#include "fastnlotk/fastNLOTools.h"

using namespace std;
using namespace fastNLO;


//________________________________________________________________________________________________________________ //
fastNLOCoeffData::fastNLOCoeffData(int NObsBin)
   : fastNLOCoeffBase(NObsBin), Nuncorrel(), UncDescr(), Ncorrel(), CorDescr(), Xcenter(), Value(),
     UncorLo(), UncorHi(), CorrLo(), CorrHi(), NErrMatrix(), matrixelement() {
   SetClassName("fastNLOCoeffData");
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffData::fastNLOCoeffData(const fastNLOCoeffBase& base) : fastNLOCoeffBase(base) {
   SetClassName("fastNLOCoeffData");
}


//________________________________________________________________________________________________________________ //
bool fastNLOCoeffData::CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet) {
   if ( c->GetIDataFlag() == 1 ) return true;
   else {
      if ( !quiet )
         say::info["fastNLOCoeffData::CheckCoeffConstants"]<<"This is not a data table! IDataFlag="<<c->GetIDataFlag()<<", but must be 1."<<endl;
      return false;
   }
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffData* fastNLOCoeffData::Clone() const {
   //! Use has to take care to delete this object later
   return new fastNLOCoeffData(*this);
}


///________________________________________________________________________________________________________________ //
void fastNLOCoeffData::Read(istream& table, int ITabVersionRead){
   fastNLOCoeffBase::ReadBase(table, ITabVersionRead);
   ReadRest(table, ITabVersionRead);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffData::ReadRest(istream& table, int ITabVersionRead){
   CheckCoeffConstants(this);
   ReadCoeffData(table);
   EndReadCoeff(table, ITabVersionRead);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffData::ReadCoeffData(istream& table){
   char buffer[5257];
   CheckCoeffConstants(this);
   table >> Nuncorrel;
   UncDescr.resize(Nuncorrel);
   table.getline(buffer,5256);
   for(int i=0;i<Nuncorrel;i++){
      table.getline(buffer,5256);
      UncDescr[i] = buffer;
      //         StripWhitespace(UncDescr[i]);
   }
   table >> Ncorrel;
   CorDescr.resize(Ncorrel);
   table.getline(buffer,5256);
   for(int i=0;i<Ncorrel;i++){
      table.getline(buffer,256);
      CorDescr[i] = buffer;
      //         StripWhitespace(CorDescr[i]);
   }
   Xcenter.resize(fNObsBins);
   Value.resize(fNObsBins);
   UncorLo.resize(fNObsBins);
   UncorHi.resize(fNObsBins);
   CorrLo.resize(fNObsBins);
   CorrHi.resize(fNObsBins);
   for(int i=0;i<fNObsBins;i++){
      table >> Xcenter[i];
      table >> Value[i];
      UncorLo[i].resize(Nuncorrel);
      UncorHi[i].resize(Nuncorrel);
      for(int j=0;j<Nuncorrel;j++){
         table >> UncorLo[i][j];
         table >> UncorHi[i][j];
      }
      CorrLo[i].resize(Ncorrel);
      CorrHi[i].resize(Ncorrel);
      for(int j=0;j<Ncorrel;j++){
         table >> CorrLo[i][j];
         table >> CorrHi[i][j];
      }
   }
   table >> NErrMatrix;
   matrixelement.resize(NErrMatrix);
   for(int i=0;i<NErrMatrix;i++){
      matrixelement[i].resize((int)pow((double)fNObsBins,2));
      for(int j=0;j<(int)pow((double)fNObsBins,2);j++){
         table >> matrixelement[i][j];
      }
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffData::Write(ostream& table, int itabversion) {
   fastNLOCoeffBase::Write(table,itabversion);
   CheckCoeffConstants(this);

   //if(IDataFlag==1){
   table << Nuncorrel << sep;
   for(int i=0;i<Nuncorrel;i++){
      table << UncDescr[i] << sep;
   }
   table << Ncorrel << sep;
   for(int i=0;i<Ncorrel;i++){
      table << CorDescr[i]  << sep;
   }
   for(int i=0;i<fNObsBins;i++){
      table << Xcenter[i] << sep;
      table << Value[i] << sep;
      for(int j=0;j<Nuncorrel;j++){
         table << UncorLo[i][j] << sep;
         table << UncorHi[i][j] << sep;
      }
      for(int j=0;j<Ncorrel;j++){
         table << CorrLo[i][j] << sep;
         table << CorrHi[i][j] << sep;
      }
   }
   table << NErrMatrix << sep;
   for(int i=0;i<NErrMatrix;i++){
      for(int j=0;j<(int)pow((double)fNObsBins,2);j++){
         table << matrixelement[i][j] << sep;
      }
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffData::Print(int iprint) const {
   if ( !(iprint < 0) ) {
      fastNLOCoeffBase::Print(iprint);
      cout << fastNLO::_DSEP20C << " fastNLO Table: CoeffData " << fastNLO::_DSEP20 << endl;
   } else {
      cout << endl << fastNLO::_CSEP20C << " fastNLO Table: CoeffData " << fastNLO::_CSEP20 << endl;
   }
   printf(" # No. of uncorr. unc. (Nuncorrel)     %d\n",Nuncorrel);
   if (Nuncorrel > 0) {fastNLOTools::PrintVector(UncDescr,"Uncorr. uncertainties (UncDescr)","#");}
   printf(" # No. of corr. unc. (Ncorrel)         %d\n",Ncorrel);
   if (Ncorrel > 0) {fastNLOTools::PrintVector(CorDescr,"Corr. uncertainties (CorDescr)","#");}
   if ( std::abs(iprint) > 0 ) {
      cout << fastNLO::_SSEP20C << " Extended information (iprint > 0) " << fastNLO::_SSEP20 << endl;
      fastNLOTools::PrintVector(Xcenter,"Data bin centers (Xcenter)","#  ");
      fastNLOTools::PrintVector(Value,"Data values (Value)","#  ");
   }
   if ( std::abs(iprint) > 1 ) {
      cout << fastNLO::_SSEP20C << " Extended information (iprint > 1) " << fastNLO::_SSEP20 << endl;
      for (int i=0; i<fNObsBins; i++) {
         // Print only for first and last observable bin
         if (i==0 || i==fNObsBins-1) {
            printf(" #     Observable bin no. %d\n",i+1);
            if (Nuncorrel > 0) {
               fastNLOTools::PrintVector(UncorLo[i],"Lower uncorr. uncertainties (UncorLo)","#    ");
               fastNLOTools::PrintVector(UncorHi[i],"Upper uncorr. uncertainties (UncorHi)","#    ");
            }
            if (Ncorrel > 0) {
               fastNLOTools::PrintVector(CorrLo[i],"Lower corr. uncertainties (CorrLo)","#    ");
               fastNLOTools::PrintVector(CorrHi[i],"Upper corr. uncertainties (CorrHi)","#    ");
            }
         }
      }
   }
   cout << fastNLO::_CSEPSC << endl;
}


//________________________________________________________________________________________________________________ //

// Erase observable bin
void fastNLOCoeffData::EraseBin(unsigned int iObsIdx) {
   debug["fastNLOCoeffData::EraseBin"]<<"Erasing table entries in CoeffData for bin index " << iObsIdx << endl;
   if ( Value.size() == 0 ) {
      say::error["EraseBin"]<<"All data bins deleted already. Aborted!" << endl;
      exit(1);
   }
   if ( Xcenter.size() != 0 ) Xcenter.erase(Xcenter.begin()+iObsIdx);
   if ( Value.size() != 0 ) Value.erase(Value.begin()+iObsIdx);
   if ( UncorLo.size() != 0 ) UncorLo.erase(UncorLo.begin()+iObsIdx);
   if ( UncorHi.size() != 0 ) UncorHi.erase(UncorHi.begin()+iObsIdx);
   if ( CorrLo.size() != 0 ) CorrLo.erase(CorrLo.begin()+iObsIdx);
   if ( CorrHi.size() != 0 ) CorrHi.erase(CorrHi.begin()+iObsIdx);
   fastNLOCoeffBase::EraseBin(iObsIdx);
}

// Catenate observable bin
void fastNLOCoeffData::CatBin(const fastNLOCoeffData& other, unsigned int iObsIdx) {
   debug["fastNLOCoeffData::CatBin"]<<"Catenating observable bin in CoeffData corresponding to bin index " << iObsIdx << endl;
   if ( Value.size() == 0 ) {
      say::error["CatBin"]<<"Initial data table is empty. Aborted!" << endl;
      exit(1);
   }
   unsigned int nold = Value.size();
   if ( Xcenter.size() != 0 ) {
      Xcenter.resize(nold+1);
      Xcenter[nold] = other.Xcenter[iObsIdx];
   }
   if ( Value.size() != 0 ) {
      Value.resize(nold+1);
      Value[nold] = other.Value[iObsIdx];
   }
   if ( UncorLo.size() != 0 ) {
      UncorLo.resize(nold+1);
      UncorLo[nold] = other.UncorLo[iObsIdx];
   }
   if ( UncorHi.size() != 0 ) {
      UncorHi.resize(nold+1);
      UncorHi[nold] = other.UncorHi[iObsIdx];
   }
   if ( CorrLo.size() != 0 ) {
      CorrLo.resize(nold+1);
      CorrLo[nold] = other.CorrLo[iObsIdx];
   }
   if ( CorrHi.size() != 0 ) {
      CorrHi.resize(nold+1);
      CorrHi[nold] = other.CorrHi[iObsIdx];
   }
   fastNLOCoeffBase::CatBin(other, iObsIdx);
}

// Multiply observable bin
void fastNLOCoeffData::MultiplyBin(unsigned int iObsIdx, double fact) {
   debug["fastNLOCoeffData::MultiplyBin"]<<"Multiplying table entries in CoeffData for bin index " << iObsIdx << " by factor " << fact << endl;
   Value[iObsIdx] *= fact;
   fastNLOCoeffBase::MultiplyBin(iObsIdx, fact);
}

//________________________________________________________________________________________________________________ //
bool fastNLOCoeffData::IsCatenable(const fastNLOCoeffData& other) const {
   //! Check for compatibility of catenating observable bins
   if ( ! ((fastNLOCoeffBase*)this)->IsCatenable(other)) return false;
   if( Nuncorrel != other.GetNuncorrel() ){
      debug["IsCatenable"]<<"Nuncorrel != other.GetNuncorrel(). Skipped."<<endl;
      return false;
   }
   if( Ncorrel != other.GetNcorrel() ){
      debug["IsCatenable"]<<"Ncorrel != other.GetNcorrel(). Skipped."<<endl;
      return false;
   }
   if( NErrMatrix != other.GetNErrMatrix() ){
      debug["IsCatenable"]<<"NErrMatrix != other.GetNErrMatrix(). Skipped."<<endl;
      return false;
   }
   info["IsCatenable"]<<"Data contributions are catenable"<<endl;
   return true;
}
