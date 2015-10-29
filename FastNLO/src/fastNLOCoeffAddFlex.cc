#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOTools.h"
#include "fastnlotk/fastNLOCoeffAddFlex.h"

using namespace std;


//________________________________________________________________________________________________________________ //
bool fastNLOCoeffAddFlex::CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet ) {
   bool ret = fastNLOCoeffAddBase::CheckCoeffConstants(c,quiet);
   if ( ret &&  c->GetNScaleDep() >= 3) return true;
   else if ( c->GetNScaleDep() < 3 ) {
      if ( !quiet )
         say::error["CheckCoeffConstants"]<<"This is not a flexible scale table. NScaleDep must be >= 3 but is NScaleDep="
                                          <<c->GetNScaleDep()<<endl;
      return false;
   }
   else return false;
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffAddFlex::fastNLOCoeffAddFlex(){
   SetClassName("fastNLOCoeffAddFlex");
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffAddFlex::fastNLOCoeffAddFlex(int NObsBin, int iLOord) : fastNLOCoeffAddBase(NObsBin){
   SetClassName("fastNLOCoeffAddFlex");
   fILOord = iLOord; // only necessary for fixing NScaleDep 3 -> 4,5
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffAddFlex::fastNLOCoeffAddFlex(const fastNLOCoeffBase& base , int iLOord ) : fastNLOCoeffAddBase(base)  {
   SetClassName("fastNLOCoeffAddFlex");
   fILOord = iLOord;
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffBase* fastNLOCoeffAddFlex::Clone() const {
   //! User has to take care to delete this object later
   return static_cast<fastNLOCoeffBase*>(new fastNLOCoeffAddFlex(*this));
}


///________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::Read(istream& table){
   fastNLOCoeffBase::ReadBase(table);
   ReadRest(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::ReadRest(istream& table){
   CheckCoeffConstants(this);
   fastNLOCoeffAddBase::ReadCoeffAddBase(table);
   ReadCoeffAddFlex(table);
   EndReadCoeff(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::ReadCoeffAddFlex(istream& table){
   CheckCoeffConstants(this);

   //  ---- order of reading... ---- //
   //    - nscalenode q2
   //    - scalenode Q
   //    - nscalenode pt
   //    - scalenode pt
   //    - simgatilde mu indep
   //    - simgatilde mu_f dep
   //    - simgatilde mu_r dep
   //    - sigmarefmixed
   //    - sigmaref scale 1
   //    - sigmaref scale 2
   // ------------------------------ //
   int nn3 = 0;

   nn3 += fastNLOTools::ReadFlexibleVector  ( ScaleNode1 , table );
   nn3 += fastNLOTools::ReadFlexibleVector  ( ScaleNode2 , table );
   //NscalenodeScale1 = ScaleNode1[0].size();
   //NscalenodeScale2 = ScaleNode2[0].size();

   nn3 += fastNLOTools::ReadFlexibleVector  ( SigmaTildeMuIndep , table , NSubproc , Nevt );
   //if ( NScaleDep==3 || fScen->ILOord!=Npow || NScaleDep==5 ){
   if ( NScaleDep==3 || NScaleDep>=5 ){
      nn3 += fastNLOTools::ReadFlexibleVector  ( SigmaTildeMuFDep , table , NSubproc , Nevt );
      nn3 += fastNLOTools::ReadFlexibleVector  ( SigmaTildeMuRDep , table , NSubproc , Nevt );
      if ( NScaleDep>=6 ){
         nn3 += fastNLOTools::ReadFlexibleVector  ( SigmaTildeMuRRDep , table , NSubproc , Nevt );
      }
      if ( NScaleDep>=7 ){
         nn3 += fastNLOTools::ReadFlexibleVector  ( SigmaTildeMuFFDep , table , NSubproc , Nevt );
         nn3 += fastNLOTools::ReadFlexibleVector  ( SigmaTildeMuRFDep , table , NSubproc , Nevt );
      }
   }
   // fixing old convention
   if ( NScaleDep == 3 ) {
      info["ReadCoeffAddFlex"]<<"This is a table with a deprecated convention (NScaleDep=3). Fixing it."<<endl;
      if (Npow!=fILOord) NScaleDep = 5;
      else NScaleDep = 3;
   }
   nn3 += fastNLOTools::ReadFlexibleVector  ( SigmaRefMixed , table , NSubproc , Nevt );
   nn3 += fastNLOTools::ReadFlexibleVector  ( SigmaRef_s1 , table , NSubproc , Nevt );
   nn3 += fastNLOTools::ReadFlexibleVector  ( SigmaRef_s2 , table , NSubproc , Nevt );
   debug["ReadCoeffAddFlex"]<<"Read "<<nn3<<" lines of flexible-scale coefficients."<<endl;

   // init table for evaluation
   fastNLOTools::ResizeFlexibleVector( PdfLcMuVar , SigmaTildeMuIndep );
   AlphasTwoPi.resize(ScaleNode1.size());
   for (unsigned int i=0; i<AlphasTwoPi.size() ; i++) {
      AlphasTwoPi[i].resize(ScaleNode1[i].size());
      for (unsigned int j=0; j<AlphasTwoPi[i].size() ; j++) {
         AlphasTwoPi[i][j].resize(ScaleNode2[i].size());
      }
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::Write(ostream& table) {
   CheckCoeffConstants(this);
   // update to latest version
   if ( NScaleDep==3 ) {
      if ( Npow==fILOord) {
         debug["Write"]<<" * Increase NScaleDep from 3 to 4, because LO!"<<endl;
         NScaleDep=4;
      }
      else if ( Npow==fILOord+1 ) {
         NScaleDep=5;
         if ( !fastNLOTools::IsEmptyVector(SigmaTildeMuRRDep) ) {
            NScaleDep=6;
            debug["Write"]<<" * Increase NScaleDep from 3 to 6 because NLO and log^2(mur) terms!"<<endl;
         }
         else
            debug["Write"]<<" * Increase NScaleDep from 3 to 5 because NLO!"<<endl;

      }
      else if ( Npow==fILOord+2 ) {
         debug["Write"]<<" * Increase NScaleDep from 3 to 6 because NNLO!"<<endl;
          NScaleDep=7;
      }
   }
   fastNLOCoeffAddBase::Write(table);

   int nn3 = 0;
   nn3 += fastNLOTools::WriteFlexibleVector( ScaleNode1 , table );
   nn3 += fastNLOTools::WriteFlexibleVector( ScaleNode2 , table );

   nn3 += fastNLOTools::WriteFlexibleVector( SigmaTildeMuIndep, table , NSubproc , Nevt);
   if ( NScaleDep==3 || NScaleDep>=5) {
      //cout<<"Write NLO FlexTable. NScaleDep="<<NScaleDep<<"\tNpow="<<Npow<<"\tfScen->ILOord="<<fScen->ILOord<<endl;
      nn3 += fastNLOTools::WriteFlexibleVector( SigmaTildeMuFDep , table , NSubproc, Nevt);
      nn3 += fastNLOTools::WriteFlexibleVector( SigmaTildeMuRDep , table , NSubproc, Nevt);
      if ( NScaleDep>=6) {
         nn3 += fastNLOTools::WriteFlexibleVector( SigmaTildeMuRRDep , table , NSubproc, Nevt);
      }
      if ( NScaleDep>=7) {
         nn3 += fastNLOTools::WriteFlexibleVector( SigmaTildeMuFFDep , table , NSubproc, Nevt);
         nn3 += fastNLOTools::WriteFlexibleVector( SigmaTildeMuRFDep , table , NSubproc, Nevt);
      }
   }
   if ( SigmaRefMixed.empty() ) fastNLOTools::ResizeVector(SigmaRefMixed,fNObsBins,NSubproc);
   if ( SigmaRef_s1.empty() )   fastNLOTools::ResizeVector(SigmaRef_s1,fNObsBins,NSubproc);
   if ( SigmaRef_s2.empty() )   fastNLOTools::ResizeVector(SigmaRef_s2,fNObsBins,NSubproc);
   nn3 += fastNLOTools::WriteFlexibleVector( SigmaRefMixed      , table , NSubproc, Nevt);
   nn3 += fastNLOTools::WriteFlexibleVector( SigmaRef_s1        , table , NSubproc, Nevt);
   nn3 += fastNLOTools::WriteFlexibleVector( SigmaRef_s2        , table , NSubproc, Nevt);

   /*
   nn3 += WriteFlexibleTable( &SigmaTildeMuIndep, table , (bool)(option & DividebyNevt) , Nevt , true );

   //if ( NScaleDep==3 || Npow!=fScen->ILOord || NScaleDep==5) {
   if ( NScaleDep==3 || NScaleDep>=5) {
      //cout<<"Write NLO FlexTable. NScaleDep="<<NScaleDep<<"\tNpow="<<Npow<<"\tfScen->ILOord="<<fScen->ILOord<<endl;
      nn3 += WriteFlexibleTable( &SigmaTildeMuFDep , table , (bool)(option & DividebyNevt) , Nevt , true );
      nn3 += WriteFlexibleTable( &SigmaTildeMuRDep , table , (bool)(option & DividebyNevt) , Nevt , true );
      if ( NScaleDep>=6) {
         nn3 += WriteFlexibleTable( &SigmaTildeMuRRDep , table , (bool)(option & DividebyNevt) , Nevt , true );
         nn3 += WriteFlexibleTable( &SigmaTildeMuFFDep , table , (bool)(option & DividebyNevt) , Nevt , true );
         nn3 += WriteFlexibleTable( &SigmaTildeMuRFDep , table , (bool)(option & DividebyNevt) , Nevt , true );
      }
   }
   if ( SigmaRefMixed.empty() ) fastNLOCoeffBase::ResizeTable(&SigmaRefMixed,fNObsBins,NSubproc);
   if ( SigmaRef_s1.empty() )   fastNLOCoeffBase::ResizeTable(&SigmaRef_s1,fNObsBins,NSubproc);
   if ( SigmaRef_s2.empty() )   fastNLOCoeffBase::ResizeTable(&SigmaRef_s2,fNObsBins,NSubproc);
   nn3 += WriteFlexibleTable( &SigmaRefMixed    , table , (bool)(option & DividebyNevt) , Nevt , true );
   nn3 += WriteFlexibleTable( &SigmaRef_s1      , table , (bool)(option & DividebyNevt) , Nevt , true );
   nn3 += WriteFlexibleTable( &SigmaRef_s2      , table , (bool)(option & DividebyNevt) , Nevt , true );
   */
   //printf("  *  fastNLOCoeffAddFlex::Write(). Wrote %d lines of v2.1 Tables.\n",nn3);
   debug["Write"]<<"Wrote "<<nn3<<" lines of v2.1 Tables."<<endl;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::Add(const fastNLOCoeffAddBase& other){
   bool ok = CheckCoeffConstants(this);
   if ( !ok ) {
      error["Add"]<<"Incompatible table."<<endl;
   }
   const fastNLOCoeffAddFlex& othflex = (const fastNLOCoeffAddFlex&) other;
   Nevt += othflex.Nevt;
   fastNLOTools::AddVectors( SigmaTildeMuIndep , othflex.SigmaTildeMuIndep );
   if ( NScaleDep==3 || NScaleDep>=5 ) {
      fastNLOTools::AddVectors( SigmaTildeMuFDep , othflex.SigmaTildeMuFDep );
      fastNLOTools::AddVectors( SigmaTildeMuRDep , othflex.SigmaTildeMuRDep );
      if (( NScaleDep>=6 || !SigmaTildeMuRRDep.empty())  // both tables contain log^2 contributions (default case)
          && (othflex.NScaleDep>=6 || !othflex.SigmaTildeMuRRDep.empty()) ) {
         fastNLOTools::AddVectors( SigmaTildeMuRRDep , othflex.SigmaTildeMuRRDep );
      }
      else if ( NScaleDep==6 && othflex.NScaleDep==5 ) { // this tables contains log^2 contributions, but the other does not
         // nothing todo.
      }
      else if ( NScaleDep==5 && othflex.NScaleDep==6 ) { // this tables does not contain log^2 contributions, but the other does !
         SigmaTildeMuRRDep = othflex.SigmaTildeMuRRDep;
         NScaleDep = 6;
      }
      if ( NScaleDep>=7 || !SigmaTildeMuFFDep.empty()) {
         fastNLOTools::AddVectors( SigmaTildeMuFFDep , othflex.SigmaTildeMuFFDep );
         fastNLOTools::AddVectors( SigmaTildeMuRFDep , othflex.SigmaTildeMuRFDep );
      }
   }
   fastNLOTools::AddVectors( SigmaRefMixed , othflex.SigmaRefMixed );
   fastNLOTools::AddVectors( SigmaRef_s1 , othflex.SigmaRef_s1 );
   fastNLOTools::AddVectors( SigmaRef_s2 , othflex.SigmaRef_s2 );
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::Clear() {
   //! Set all elements of sigma tilde to zero
   fastNLOCoeffAddBase::Clear();
   fastNLOTools::ClearVector(SigmaTildeMuIndep);
   fastNLOTools::ClearVector(SigmaTildeMuFDep);
   fastNLOTools::ClearVector(SigmaTildeMuRDep);
   fastNLOTools::ClearVector(SigmaTildeMuRRDep);
   fastNLOTools::ClearVector(SigmaTildeMuFFDep);
   fastNLOTools::ClearVector(SigmaTildeMuRFDep);
   fastNLOTools::ClearVector(SigmaRefMixed);
   fastNLOTools::ClearVector(SigmaRef_s1);
   fastNLOTools::ClearVector(SigmaRef_s2);
}

//________________________________________________________________________________________________________________ // 
bool  fastNLOCoeffAddFlex::IsCompatible(const fastNLOCoeffAddFlex& other) const {
   //! Check for compatibility for merging/adding of two contributions 
   if ( ! ((fastNLOCoeffAddBase*)this)->IsCompatible(other)) return false;
   for ( int i=0 ; i<fNObsBins ; i++ ){ 
      if ( GetNScaleNode1(i) != other.GetNScaleNode1(i) ) {
	 say::warn["fastNLOCoeffAddFlex::IsCompatible"]<<"Incompatible number of scale nodes found."<<endl;
	 return false;
      }
      if ( GetNScaleNode2(i) != other.GetNScaleNode2(i) ) {
	 say::warn["fastNLOCoeffAddFlex::IsCompatible"]<<"Incompatible number of scale nodes found."<<endl;
	 return false;
      }
      for ( unsigned int is1 = 0 ; is1<GetNScaleNode1(i) ; is1++ ) {
	 if ( GetScaleNode1(i,is1) != other.GetScaleNode1(i,is1) ) {
	    say::warn["fastNLOCoeffAddFlex::IsCompatible"]<<"Incompatible scale1 node found."<<endl;
	    return false;
	 }
      }
      for ( unsigned int is2 = 0 ; is2<GetNScaleNode2(i) ; is2++ ) {
	 if ( GetScaleNode2(i,is2) != other.GetScaleNode2(i,is2) ) {
	    say::warn["fastNLOCoeffAddFlex::IsCompatible"]<<"Incompatible scale2 node found."<<endl;
	    return false;
	 }
      }
   }
   return true;
}

//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::NormalizeCoefficients(){
   //!< Set number of events to 1 and normalize coefficients accordingly.
   //! This means, that the information about the
   //! number of events is essentially lost
   MultiplyCoefficientsByConstant(1./Nevt);
   Nevt = 1;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::MultiplyCoefficientsByConstant(double coef) {
   const bool WithMuR = !SigmaTildeMuFDep.empty();
   const bool WithMuRR = !SigmaTildeMuRRDep.empty();
   const bool WithMuFF = !SigmaTildeMuFFDep.empty();
   for (unsigned int i=0; i<SigmaTildeMuIndep.size(); i++) {
      int nxmax = GetNxmax(i);
      for (unsigned int jS1=0; jS1<GetNScaleNode1(i); jS1++) {
         for (unsigned int kS2=0; kS2<GetNScaleNode2(i); kS2++) {
            for (int x=0; x<nxmax; x++) {
               for (int n=0; n<GetNSubproc(); n++) {
                  SigmaTildeMuIndep[i][x][jS1][kS2][n] *= coef;
                  //if ( GetNScaleDep() >= 5 ) {
                  if (WithMuR) {
                     SigmaTildeMuFDep [i][x][jS1][kS2][n] *= coef;
                     SigmaTildeMuRDep [i][x][jS1][kS2][n] *= coef;
                     //if ( GetNScaleDep() >= 6 ) {
                     if (WithMuRR) {
                        SigmaTildeMuRRDep [i][x][jS1][kS2][n] *= coef;
                     }
                     if (WithMuFF) {
                        SigmaTildeMuFFDep [i][x][jS1][kS2][n] *= coef;
                        SigmaTildeMuRFDep [i][x][jS1][kS2][n] *= coef;
                     }
                  }
               }
            }
         }
      }
   }
}

//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::Print() const {
   fastNLOCoeffAddBase::Print();
   printf(" **************** FastNLO Table: fastNLOCoeffAddFlex ****************\n");
   printf(" B   NscalenodeScale1              %lu\n",ScaleNode1[0].size());
   printf(" B   NscalenodeScale2              %lu\n",ScaleNode2[0].size());
   printf(" *******************************************************\n");
}


//________________________________________________________________________________________________________________ //
