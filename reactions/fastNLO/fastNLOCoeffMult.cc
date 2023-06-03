#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOCoeffMult.h"
#include "fastnlotk/fastNLOTools.h"

using namespace std;
using namespace fastNLO;


//________________________________________________________________________________________________________________ //
fastNLOCoeffMult::fastNLOCoeffMult(int NObsBin)
   : fastNLOCoeffBase(NObsBin), Nuncorrel(), UncDescr(), Ncorrel(), CorDescr(),
     UncorLo(), UncorHi(), CorrLo(), CorrHi(), fact() {
   SetClassName("fastNLOCoeffMult");
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffMult::fastNLOCoeffMult(const fastNLOCoeffBase& base) : fastNLOCoeffBase(base) {
   SetClassName("fastNLOCoeffMult");
}


//________________________________________________________________________________________________________________ //
bool fastNLOCoeffMult::CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet)  {
   if ( c->GetIAddMultFlag()==0 && c->GetIDataFlag()==0 ) {
      // Additive contribution
      return false;
   } else if ( c->GetIAddMultFlag()==1 && c->GetIDataFlag()==0 ) {
      // Multiplicative contribution
      return true;
   } else if ( c->GetIAddMultFlag()==0 && c->GetIDataFlag()==1 ) {
      // Data contribution
      return false;
   } else {
      // Unknown contribution
      say::error["fastNLOCoeffMult::CheckCoeffConstants"]
         << "Unknown contribution type, aborting! "
         << "IAddMultFlag = " << c->GetIAddMultFlag()
         << ", IDataFlag ="   << c->GetIDataFlag() <<endl;
      exit(1);
   }
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffMult* fastNLOCoeffMult::Clone() const {
   //! Use has to take care to delete this object later
   return new fastNLOCoeffMult(*this);
}


///________________________________________________________________________________________________________________ //
void fastNLOCoeffMult::Read(istream& table, int ITabVersionRead){
   fastNLOCoeffBase::ReadBase(table, ITabVersionRead);
   ReadRest(table, ITabVersionRead);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffMult::ReadRest(istream& table, int ITabVersionRead){
   CheckCoeffConstants(this);
   ReadCoeffMult(table);
   EndReadCoeff(table, ITabVersionRead);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffMult::ReadCoeffMult(istream& table){
   char buffer[5257];
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
      table.getline(buffer,5256);
      CorDescr[i] = buffer;
      //         StripWhitespace(CorDescr[i]);
   }
   fact.resize(fNObsBins);
   UncorLo.resize(fNObsBins);
   UncorHi.resize(fNObsBins);
   CorrLo.resize(fNObsBins);
   CorrHi.resize(fNObsBins);
   for(int i=0;i<fNObsBins;i++){
      table >> fact[i];
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
   // end of IAddMultFlag==1
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffMult::Write(ostream& table, int itabversion) {
   fastNLOCoeffBase::Write(table,itabversion);
   CheckCoeffConstants(this);
   table << Nuncorrel << sep;
   for(int i=0;i<Nuncorrel;i++){
      table << UncDescr[i]  << sep;
   }
   table << Ncorrel << sep;
   for(int i=0;i<Ncorrel;i++){
      table << CorDescr[i]  << sep;
   }
   for(int i=0;i<fNObsBins;i++){
      table << fact[i] << sep;
      for(int j=0;j<Nuncorrel;j++){
         table << UncorLo[i][j] << sep;
         table << UncorHi[i][j] << sep;
      }
      for(int j=0;j<Ncorrel;j++){
         table << CorrLo[i][j] << sep;
         table << CorrHi[i][j] << sep;
      }
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffMult::Print(int iprint) const {
   if ( !(iprint < 0) ) {
      fastNLOCoeffBase::Print(iprint);
      cout << fastNLO::_DSEP20C << " fastNLO Table: CoeffMult " << fastNLO::_DSEP20 << endl;
   } else {
      cout << endl << fastNLO::_CSEP20C << " fastNLO Table: CoeffMult " << fastNLO::_CSEP20 << endl;
   }
   double minfact = *min_element(fact.begin(),fact.end());
   double maxfact = *max_element(fact.begin(),fact.end());
   printf(" # Minimal correction factor (fact[])  %f\n",minfact);
   printf(" # Maximal correction factor (fact[])  %f\n",maxfact);
   printf(" # No. of uncorr. unc. (Nuncorrel)     %d\n",Nuncorrel);
   if (Nuncorrel > 0) {fastNLOTools::PrintVector(UncDescr,"Uncorr. uncertainties (UncDescr)","#");}
   printf(" # No. of corr. unc. (Ncorrel)         %d\n",Ncorrel);
   if (Ncorrel > 0) {fastNLOTools::PrintVector(CorDescr,"Corr. uncertainties (CorDescr)","#");}
   if ( std::abs(iprint) > 1 ) {
      cout << fastNLO::_SSEP20C << " Extended information (iprint > 1) " << fastNLO::_SSEP20 << endl;
      fastNLOTools::PrintVector(fact,"Correction factors (fact)","#    ");
   }
   if ( std::abs(iprint) > 2 ) {
      cout << fastNLO::_SSEP20C << " Extended information (iprint > 2) " << fastNLO::_SSEP20 << endl;
      for (int i=0; i<fNObsBins; i++) {
         // Print only for first and last observable bin
         if (i==0 || i==fNObsBins-1) {
            printf(" #       Observable bin no. %d\n",i+1);
            if (Nuncorrel > 0) {
               fastNLOTools::PrintVector(UncorLo[i],"Lower uncorr. uncertainties (UncorLo)","#      ");
               fastNLOTools::PrintVector(UncorHi[i],"Upper uncorr. uncertainties (UncorHi)","#      ");
            }
            if (Ncorrel > 0) {
               fastNLOTools::PrintVector(CorrLo[i],"Lower corr. uncertainties (CorrLo)","#      ");
               fastNLOTools::PrintVector(CorrHi[i],"Upper corr. uncertainties (CorrHi)","#      ");
            }
         }
      }
   }
   cout << fastNLO::_CSEPSC << endl;
}


//________________________________________________________________________________________________________________ //
// Erase observable bin
void fastNLOCoeffMult::EraseBin(unsigned int iObsIdx) {
   debug["fastNLOCoeffMult::EraseBin"]<<"Erasing table entries in CoeffMult for bin index " << iObsIdx << endl;
   if ( fact.size() == 0 ) {
      say::error["EraseBin"]<<"All multiplicative bins deleted already. Aborted!" << endl;
      exit(1);
   }
   if ( fact.size() != 0 ) fact.erase(fact.begin()+iObsIdx);
   if ( UncorLo.size() != 0 ) UncorLo.erase(UncorLo.begin()+iObsIdx);
   if ( UncorHi.size() != 0 ) UncorHi.erase(UncorHi.begin()+iObsIdx);
   if ( CorrLo.size() != 0 ) CorrLo.erase(CorrLo.begin()+iObsIdx);
   if ( CorrHi.size() != 0 ) CorrHi.erase(CorrHi.begin()+iObsIdx);
   fastNLOCoeffBase::EraseBin(iObsIdx);
}

// Catenate observable bin
void fastNLOCoeffMult::CatBin(const fastNLOCoeffMult& other, unsigned int iObsIdx) {
   debug["fastNLOCoeffMult::CatBin"]<<"Catenating observable bin in CoeffMult corresponding to bin index " << iObsIdx << endl;
   if ( fact.size() == 0 ) {
      say::error["CatBin"]<<"Initial multiplicative table is empty. Aborted!" << endl;
      exit(1);
   }
   unsigned int nold = fact.size();
   if ( fact.size() != 0 ) {
      fact.resize(nold+1);
      fact[nold] = other.fact[iObsIdx];
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
void fastNLOCoeffMult::MultiplyBin(unsigned int iObsIdx, double nfact) {
   debug["fastNLOCoeffMult::MultiplyBin"]<<"Multiplying table entries. Nothing to be done in CoeffMult." << endl;
}

//________________________________________________________________________________________________________________ //
bool fastNLOCoeffMult::IsCatenable(const fastNLOCoeffMult& other) const {
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
   info["IsCatenable"]<<"Multiplicable contributions are catenable"<<endl;
   return true;
}
