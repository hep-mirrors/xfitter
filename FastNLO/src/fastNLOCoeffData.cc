#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOCoeffData.h"

using namespace std;
using namespace fastNLO;


//________________________________________________________________________________________________________________ //
bool fastNLOCoeffData::CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet) {
   if ( c->GetIDataFlag() == 1 ) return true;
   else {
      if ( !quiet )
	 say::error["fastNLOCoeffData::CheckCoeffConstants"]<<"This is not a data table! IDataFlag="<<c->GetIDataFlag()<<", but must be 1."<<endl;
      return false;
   }
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffData::fastNLOCoeffData(){
   SetClassName("fastNLOCoeffData");
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffData::fastNLOCoeffData(int NObsBin) : fastNLOCoeffBase(NObsBin) {
   SetClassName("fastNLOCoeffData");
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffData::fastNLOCoeffData(const fastNLOCoeffBase& base) : fastNLOCoeffBase(base) {
   SetClassName("fastNLOCoeffData");
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffBase* fastNLOCoeffData::Clone() const {
   //! Use has to take care to delete this object later
   return static_cast<fastNLOCoeffBase*>(new fastNLOCoeffData(*this));
}


///________________________________________________________________________________________________________________ //
void fastNLOCoeffData::Read(istream& table){
   fastNLOCoeffBase::ReadBase(table);
   ReadRest(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffData::ReadRest(istream& table){
   CheckCoeffConstants(this);
   ReadCoeffData(table);
   EndReadCoeff(table);
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
void fastNLOCoeffData::Write(ostream& table) {
   fastNLOCoeffBase::Write(table);
   CheckCoeffConstants(this);

   //if(IDataFlag==1){
   table << Nuncorrel << endl;
   for(int i=0;i<Nuncorrel;i++){
      table << UncDescr[i] << endl;
   }
   table << Ncorrel << endl;
   for(int i=0;i<Ncorrel;i++){
      table << CorDescr[i]  << endl;
   }
   for(int i=0;i<fNObsBins;i++){
      table << Xcenter[i] << endl;
      table << Value[i] << endl;
      for(int j=0;j<Nuncorrel;j++){
	 table << UncorLo[i][j] << endl;
	 table << UncorHi[i][j] << endl;
      }
      for(int j=0;j<Ncorrel;j++){
	 table << CorrLo[i][j] << endl;
	 table << CorrHi[i][j] << endl;
      }
   }
   table << NErrMatrix << endl;
   for(int i=0;i<NErrMatrix;i++){
      for(int j=0;j<(int)pow((double)fNObsBins,2);j++){
	 table << matrixelement[i][j] << endl;
      }
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffData::Print() const {
   fastNLOCoeffBase::Print();
   printf("\n **************** FastNLO Table: CoeffData ****************\n\n");
   printf(" B   Nuncorrel                     %d\n",Nuncorrel);
   printf(" B   Ncorrel                       %d\n",Ncorrel);
   printf(" B   NErrMatrix                    %d\n",NErrMatrix);
   printf(" B   fastNLOCoeffData::Print(). some more output could be printed here (IDataFlag==1).\n");
   printf("\n *******************************************************\n\n");
}


//________________________________________________________________________________________________________________ //
