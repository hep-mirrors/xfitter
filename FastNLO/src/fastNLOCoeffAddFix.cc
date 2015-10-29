#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOTools.h"
#include "fastnlotk/fastNLOCoeffAddFix.h"
#include "fastnlotk/speaker.h"

using namespace std;


//________________________________________________________________________________________________________________ //
bool fastNLOCoeffAddFix::CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet) {
   bool ret = fastNLOCoeffAddBase::CheckCoeffConstants(c,quiet);
   if ( ret && c->GetNScaleDep() == 0 ) return true;
   else if ( c->GetNScaleDep() >= 3 ) {
      if ( !quiet)
         say::error["fastNLOCoeffAddFix::CheckCoeffConstants"]
            <<"This is not a fixed order v2.0  table. NScaleDep must be equal 0 but is NScaleDep="
            <<c->GetNScaleDep()<<endl;
      return false;
   }
   else return false;
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffAddFix::fastNLOCoeffAddFix(){
   SetClassName("fastNLOCoeffAddFix");
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffAddFix::fastNLOCoeffAddFix(int NObsBin) : fastNLOCoeffAddBase(NObsBin) {
   SetClassName("fastNLOCoeffAddFix");
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffAddFix::fastNLOCoeffAddFix(const fastNLOCoeffBase& base) : fastNLOCoeffAddBase(base) {
   SetClassName("fastNLOCoeffAddFix");
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffBase* fastNLOCoeffAddFix::Clone() const {
   //! Use has to take care to delete this object later
   return static_cast<fastNLOCoeffBase*>(new fastNLOCoeffAddFix(*this));
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Read(istream& table){
   fastNLOCoeffBase::ReadBase(table);
   ReadRest(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::ReadRest(istream& table){
   CheckCoeffConstants(this);
   fastNLOCoeffAddBase::ReadCoeffAddBase(table);
   ReadCoeffAddFix(table);
   EndReadCoeff(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::ReadCoeffAddFix(istream& table){
   CheckCoeffConstants(this);

   Nscalevar.resize(NScaleDim);
   vector<int> Nscalenode(NScaleDim);
   for(int i=0;i<NScaleDim;i++){
      table >> Nscalevar[i];
      table >> Nscalenode[i];
   }
   //    printf("  *  fastNLOCoeffAddFix::Read().bins %d, NScalevar[0] %d, Nscalenode[0] %d,  NScaleDim %d  \n",
   //    fNObsBins, Nscalevar[0] , Nscalenode[0] , NScaleDim );
   ScaleFac.resize(NScaleDim);
   for(int i=0;i<NScaleDim;i++){
      ScaleFac[i].resize(Nscalevar[i]);
      for(int j=0;j<Nscalevar[i];j++){
         table >> ScaleFac[i][j];
      }
   }
   //printf("  *  fastNLOCoeffAddFix::Read().bins %d, NScalevar[0] %d, Nscalenode[0] %d, ScaleFac[0][0] %d,  NScaleDim %d  \n",
   //fNObsBins, Nscalevar[0] , Nscalenode[0] , ScaleFac[0][0], NScaleDim );
   fastNLOTools::ResizeVector( ScaleNode , fNObsBins, 1 , Nscalevar[0] , Nscalenode[0] ); // should work, since NScaleDim==1, but is not yet tested for 100%
   int nsn = fastNLOTools::ReadVector( ScaleNode , table );
   debug["fastNLOCoeffAddFix::Read()"]<<"Read "<<nsn<<" lines of ScaleNode."<<endl;

   ResizeSigmaTilde();
   ResizePdfLC();
   ResizePdfSplLC();
   int nst = fastNLOTools::ReadVector( SigmaTilde , table , Nevt);
   debug["fastNLOCoeffAddFix::Read()"]<<"Read "<<nsn+nst<<" lines of fastNLO v2 tables."<<endl;

   // prepare members for evaluation
   fastNLOTools::ResizeVector(AlphasTwoPi_v20 , fNObsBins, GetTotalScalenodes());
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::ResizeSigmaTilde(){
   //! resize SigmaTilde and ensure that all entries are empty
   SigmaTilde.resize(fNObsBins);
   for( int i=0 ; i<fNObsBins ; i++ ){
      int nxmax = GetNxmax(i);
      SigmaTilde[i].resize(GetTotalScalevars());
      for( int k=0 ; k<GetTotalScalevars() ; k++ ){
         SigmaTilde[i][k].resize(GetTotalScalenodes());
         for( int l=0 ; l<GetTotalScalenodes() ; l++ ){
            //ResizeVector(SigmaTilde[i][k][l],nxmax,NSubproc);
            SigmaTilde[i][k][l].resize(nxmax);
            for( int m=0 ; m<nxmax ; m++ ){
               SigmaTilde[i][k][l][m].resize(NSubproc);
               for( int n=0 ; n<NSubproc ; n++ ){
                  SigmaTilde[i][k][l][m][n] = 0.;
               }
            }
         }
      }
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::ResizePdfLC(){
   //! resize PdfLC
   PdfLc.resize(fNObsBins);
   for( int i=0 ; i<fNObsBins ; i++ ){
      int nxmax = GetNxmax(i);
      int totalscalenodes = GetTotalScalenodes();
      PdfLc[i].resize(totalscalenodes);
      for( int l=0 ; l<totalscalenodes ; l++ ){
         fastNLOTools::ResizeVector(PdfLc[i][l],nxmax,NSubproc);
      }
   }
}

//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::ResizePdfSplLC(){
   //! resize PdfSplLC
   PdfSplLc1.resize(fNObsBins);
   PdfSplLc2.resize(fNObsBins);
   for( int i=0 ; i<fNObsBins ; i++ ){
      int nxmax = GetNxmax(i);
      int totalscalenodes = GetTotalScalenodes();
      PdfSplLc1[i].resize(totalscalenodes);
      PdfSplLc2[i].resize(totalscalenodes);
      for( int l=0 ; l<totalscalenodes ; l++ ){
         fastNLOTools::ResizeVector(PdfSplLc1[i][l],nxmax,NSubproc);
         fastNLOTools::ResizeVector(PdfSplLc2[i][l],nxmax,NSubproc);
      }
   }
}



//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Write(ostream& table){
   //! Write coefficient table to disk (ostream)
   CheckCoeffConstants(this);
   fastNLOCoeffAddBase::Write(table);

   for(int i=0;i<NScaleDim;i++){
      table << Nscalevar[i] << endl;
      table << GetNScaleNode() << endl;
   }
   fastNLOTools::WriteVector( ScaleFac , table );
   int nsn = fastNLOTools::WriteVector( ScaleNode , table );
   int nst = fastNLOTools::WriteVector( SigmaTilde , table , Nevt);
   info["Write"]<<"Wrote "<<nst+nsn<<" lines into fastNLO table."<<endl;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Add(const fastNLOCoeffAddBase& other){
   //! Add another coefficient table to this table
   bool ok = CheckCoeffConstants(this);
   if ( !ok ) {
      error["Add"]<<"Cannot add tables."<<endl;
      return;
   }
   const fastNLOCoeffAddFix& othfix = (const fastNLOCoeffAddFix&)other;
   Nevt += othfix.Nevt;
   fastNLOTools::AddVectors( SigmaTilde , othfix.SigmaTilde);
}



//________________________________________________________________________________________________________________ //
int fastNLOCoeffAddFix::GetTotalScalevars() const {
   //! Get nuber of scale-variations
   int totalscalevars=1;
   for(int scaledim=0;scaledim<NScaleDim;scaledim++){
      totalscalevars *= Nscalevar[scaledim];
   }
   return totalscalevars;
}


//________________________________________________________________________________________________________________ //
int fastNLOCoeffAddFix::GetTotalScalenodes() const {
   //! Get number of scale nodes
   if ( !ScaleNode.empty() ) return ScaleNode[0][0][0].size();
   else return 0;
   //    int totalscalenodes=1;
   //    for(int scaledim=0;scaledim<NScaleDim;scaledim++){
   //       if ( !ScaleNode.empty() )
   //  totalscalenodes *= ScaleNode[0][0][0].size();
   //       //Nscalenode[scaledim];
   //    }
   //    return totalscalenodes;
}


//________________________________________________________________________________________________________________ //
bool  fastNLOCoeffAddFix::IsCompatible(const fastNLOCoeffAddFix& other) const {
   //! Check for compatibility for merging/adding of two contributions
   if ( ! ((fastNLOCoeffAddBase*)this)->IsCompatible(other)) return false;
   if ( GetNScaleNode() != other.GetNScaleNode() ) {
      say::warn["IsCompatible"]<<"Incompatible number of scale nodes found."<<endl;                                                                                     
      return false;
   }
   if ( GetNScalevar() != other.GetNScalevar() ) {
      say::warn["IsCompatible"]<<"Incompatible number of scale variations found."<<endl;                                                                                     
      return false;
   }
   if ( GetAvailableScaleFactors()[GetNScalevar()-1] != other.GetAvailableScaleFactors()[GetNScalevar()-1] ) {
      say::warn["IsCompatible"]<<"Incompatible scale variations found."<<endl;                                                                                     
      return false;
   }
   for ( int i=0 ; i<fNObsBins ; i++ ){
      for ( int is=0 ; is<GetNScaleNode() ; is++ ){
	 if ( GetScaleNode(i,0,is) != other.GetScaleNode(i,0,is) ) {
	    say::warn["IsCompatible"]<<"Incompatible scale node found."<<endl;                                                                                     
	    return false;
	 }
      }
   }
   return true; 
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Clear() {
   //! Set all elelments of SigmaTilde to zero.
   fastNLOCoeffAddBase::Clear();
   fastNLOTools::ClearVector(SigmaTilde);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::NormalizeCoefficients(){
   //!< Set number of events to 1 and normalize coefficients accordingly.
   //! This means, that the information about the
   //! number of events is essentially lost
   MultiplyCoefficientsByConstant(1./Nevt);
   Nevt = 1;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::MultiplyCoefficientsByConstant(double coef) {
   for (unsigned int i=0; i<SigmaTilde.size(); i++) {
      for (unsigned int s=0 ; s<SigmaTilde[i].size() ; s++) {
         for (unsigned int x=0 ; x<SigmaTilde[i][s].size() ; x++) {
            for (unsigned int l=0 ; l<SigmaTilde[i][s][x].size() ; l++) {
               for (unsigned int m=0 ; m<SigmaTilde[i][s][x][m].size() ; m++) {
                  SigmaTilde[i][s][x][l][m] *= coef;
               }
            }
         }
      }
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Print() const {
   fastNLOCoeffAddBase::Print();
   printf(" **************** FastNLO Table: fastNLOCoeffAddFix ****************\n");
   for(int i=0;i<NScaleDim;i++){
      printf(" B    - Nscalenode[%d]              %d\n",i,GetNScaleNode());
      printf(" B    - Nscalevar[%d]               %d\n",i,Nscalevar[i]);
      for(int j=0;j<Nscalevar[i];j++){
         printf(" B    -  - ScaleFac[%d][%d]          %6.4f\n",i,j,ScaleFac[i][j]);
      }
   }
   printf(" B   No printing of ScaleNode implemented yet.\n");
   printf(" B   No printing of SigmaTilde implemented yet.\n");
   printf(" *******************************************************\n");
}


//________________________________________________________________________________________________________________ //
