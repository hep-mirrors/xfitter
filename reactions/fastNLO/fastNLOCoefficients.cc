#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOCoefficients.h"

using namespace std;
using namespace fastNLO;

fastNLOCoefficients::fastNLOCoefficients(){
}

fastNLOCoefficients::fastNLOCoefficients(int NObsBin, int iLOord){
   fNObsBins = NObsBin;
   fILOord = iLOord; // only necessary for fixing NScaleDep 3 -> 4,5
   NScaleDep = 0;
}


int fastNLOCoefficients::Read(istream *table){
   table->peek();
   if (table->eof()){
      printf("fastNLOCoefficients::Read: Cannot read from file.\n");
      return(2);
   }

   int key = 0;
   *table >> key;
   if(key != tablemagicno){
      printf("fastNLOCoefficients::Read: At beginning of block found %d instead of %d.\n",key,tablemagicno);
      return 1;
   };

   *table >> IXsectUnits;
   *table >> IDataFlag;
   *table >> IAddMultFlag;
   *table >> IContrFlag1;
   *table >> IContrFlag2;
   // KR: Let's simply drop IContrFlag3, which in pp scenarios was always 0 anyway and
   // KR: reuse this table line for NScaleDep!
   // *table >> IContrFlag3;    // IContrFlag3 is written here in v2.0 and v2.1 but  not in v2.0+
   // in v2.1. IContrFlag3 will be reintroduces again, and NScaleDep will be stored later in the table
   //   IContrFlag3 = 0;
   *table >> NScaleDep;
   int NContrDescr;
   *table >> NContrDescr;
   //   printf("  *  fastNLOCoefficients::Read().  IDataFlag: %d, IAddMultFlag: %d, IContrFlag1: %d, IContrFlag2: %d,, NScaleDep: %d\n",IDataFlag,IAddMultFlag,IContrFlag1,IContrFlag2,NScaleDep );
   CtrbDescript.resize(NContrDescr);
   char buffer[257];
   table->getline(buffer,256);
   for(int i=0;i<NContrDescr;i++){
      table->getline(buffer,256);
      CtrbDescript[i] = buffer;
      //      StripWhitespace(CtrbDescript[i]);
   }
   int NCodeDescr;
   *table >> NCodeDescr;
   CodeDescript.resize(NCodeDescr);
   table->getline(buffer,256);
   for(int i=0;i<NCodeDescr;i++){
      table->getline(buffer,256);
      CodeDescript[i] = buffer;
      //      StripWhitespace(CodeDescript[i]);
   }
   if(IDataFlag==1){
      *table >> Nuncorrel;
      UncDescr.resize(Nuncorrel);
      table->getline(buffer,256);
      for(int i=0;i<Nuncorrel;i++){
         table->getline(buffer,256);
         UncDescr[i] = buffer;
         //         StripWhitespace(UncDescr[i]);
      }
      *table >> Ncorrel;
      CorDescr.resize(Ncorrel);
      table->getline(buffer,256);
      for(int i=0;i<Ncorrel;i++){
         table->getline(buffer,256);
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
         *table >> Xcenter[i];
         *table >> Value[i];
         UncorLo[i].resize(Nuncorrel);
         UncorHi[i].resize(Nuncorrel);
         for(int j=0;j<Nuncorrel;j++){
            *table >> UncorLo[i][j];
            *table >> UncorHi[i][j];
         }
         CorrLo[i].resize(Ncorrel);
         CorrHi[i].resize(Ncorrel);
         for(int j=0;j<Ncorrel;j++){
            *table >> CorrLo[i][j];
            *table >> CorrHi[i][j];
         }
      }
      *table >> NErrMatrix;
      matrixelement.resize(NErrMatrix);
      for(int i=0;i<NErrMatrix;i++){
         matrixelement[i].resize((int)pow((double)fNObsBins,2));
         for(int j=0;j<(int)pow((double)fNObsBins,2);j++){
            *table >> matrixelement[i][j];
         }
      }
   }// end of IDataFlag==1
   if(IAddMultFlag==1){
      *table >> Nuncorrel;
      UncDescr.resize(Nuncorrel);
      table->getline(buffer,256);
      for(int i=0;i<Nuncorrel;i++){
         table->getline(buffer,256);
         UncDescr[i] = buffer;
         //         StripWhitespace(UncDescr[i]);
      }
      *table >> Ncorrel;
      CorDescr.resize(Ncorrel);
      table->getline(buffer,256);
      for(int i=0;i<Ncorrel;i++){
         table->getline(buffer,256);
         CorDescr[i] = buffer;
         //         StripWhitespace(CorDescr[i]);
      }
      fact.resize(fNObsBins);
      UncorLo.resize(fNObsBins);
      UncorHi.resize(fNObsBins);
      CorrLo.resize(fNObsBins);
      CorrHi.resize(fNObsBins);
      for(int i=0;i<fNObsBins;i++){
         *table >> fact[i];
         UncorLo[i].resize(Nuncorrel);
         UncorHi[i].resize(Nuncorrel);
         for(int j=0;j<Nuncorrel;j++){
            *table >> UncorLo[i][j];
            *table >> UncorHi[i][j];
         }
         CorrLo[i].resize(Ncorrel);
         CorrHi[i].resize(Ncorrel);
         for(int j=0;j<Ncorrel;j++){
            *table >> CorrLo[i][j];
            *table >> CorrHi[i][j];
         }
      }
   }// end of IAddMultFlag==1

   if(!(IDataFlag==1) && !(IAddMultFlag==1)){
      *table >> IRef;
      *table >> IScaleDep;
      *table >> Nevt;
      *table >> Npow;
      *table >> NPDF;
      if(NPDF>0){
         NPDFPDG.resize(NPDF);
         for(int i=0;i<NPDF;i++){
            *table >>  NPDFPDG[i];
         }
      }
      *table >> NPDFDim;
      *table >> NFragFunc;
      if(NFragFunc>0){
         NFFPDG.resize(NFragFunc);
         for(int i=0;i<NFragFunc;i++){
            *table >>  NFFPDG[i];
         }
      }
      *table >> NFFDim;
      *table >> NSubproc;
      *table >> IPDFdef1;
      *table >> IPDFdef2;
      *table >> IPDFdef3;
      //printf("  *  fastNLOCoefficients::Read(). IRef : %d, IScaleDep: %d, Nevt: %d, Npow: %d, NPDF: %d, NPDFDim: %d\n", IRef ,IScaleDep  ,Nevt  , Npow ,NPDF , NPDFDim  );

      if(IPDFdef1==0){
         for(int i=0;i<NSubproc;i++){
            // Missing: linear PDF combinations for IPDFdef1=0
            if(NPDF==1){
            }else{
               if(NPDF==2){
               }
            }
         }
      }
      //Nxtot1.resize(fNObsBins);
      XNode1.resize(fNObsBins);
      for(int i=0;i<fNObsBins;i++){
         int xtot;
         *table >> xtot;
         //*table >> Nxtot1[i];
         //XNode1[i].resize(Nxtot1[i]);
         XNode1[i].resize(xtot);
         for(int j=0;j<xtot;j++){
            *table >> XNode1[i][j];
         }
      }
      if(NPDFDim==2){
         //Nxtot2.resize(fNObsBins);
         XNode2.resize(fNObsBins);
         for(int i=0;i<fNObsBins;i++){
            int xtot;
            *table >> xtot;
            XNode2[i].resize(xtot);
            //*table >> Nxtot2[i];
            //XNode2[i].resize(Nxtot2[i]);
            for(int j=0;j<xtot;j++){
               *table >> XNode2[i][j];
            }
         }
      }
      if(NFragFunc>0){
         Nztot.resize(fNObsBins);
         ZNode.resize(fNObsBins);
         for(int i=0;i<fNObsBins;i++){
            *table >> Nztot[i];
            ZNode[i].resize(Nztot[i]);
            for(int j=0;j<Nztot[i];j++){
               *table >> ZNode[i][j];
            }
         }
      }

      *table >> NScales;
      *table >> NScaleDim;
      Iscale.resize(NScales);
      for(int i=0;i<NScales;i++){
         *table >> Iscale[i];
      }
      int NscaleDescript;
      ScaleDescript.resize(NScaleDim);
      for(int i=0;i<NScaleDim;i++){
         *table >> NscaleDescript;
         ScaleDescript[i].resize(NscaleDescript);
         table->getline(buffer,256);
         for(int j=0;j<NscaleDescript;j++){
            table->getline(buffer,256);
            ScaleDescript[i][j] = buffer;
            //            StripWhitespace(ScaleDescript[i][j]);
         }
      }


      //! v2.1 store NScaleDep here.
      //! v2.1 *table >> NScaleDep;

      if ( NScaleDep < 3 ) {
         Nscalevar.resize(NScaleDim);
         Nscalenode.resize(NScaleDim);
         for(int i=0;i<NScaleDim;i++){
            *table >> Nscalevar[i];
            *table >> Nscalenode[i];
         }

         //      printf("  *  fastNLOCoefficients::Read().bins %d, NScalevar[0] %d, Nscalenode[0] %d,  NScaleDim %d  \n",
         //      fNObsBins, Nscalevar[0] , Nscalenode[0] , NScaleDim );


         ScaleFac.resize(NScaleDim);
         for(int i=0;i<NScaleDim;i++){
            ScaleFac[i].resize(Nscalevar[i]);
            for(int j=0;j<Nscalevar[i];j++){
               *table >> ScaleFac[i][j];
            }
         }

         //printf("  *  fastNLOCoefficients::Read().bins %d, NScalevar[0] %d, Nscalenode[0] %d, ScaleFac[0][0] %d,  NScaleDim %d  \n",
         //fNObsBins, Nscalevar[0] , Nscalenode[0] , ScaleFac[0][0], NScaleDim );
         ResizeTable( &ScaleNode , fNObsBins, 1 , Nscalevar[0] , Nscalenode[0] ); // should work, since NScaleDim==1, but is not yet tested for 100%
         int nsn = ReadTable  ( &ScaleNode , table );
         //printf("  *  fastNLOCoefficients::Read(). Read %d lines of ScaleNode.\n",nsn);

         int XmaxFromI[1] = {0};
         //printf(" &SigmaTilde  %i  %i  %i  *%i  %i\n",
         //fNObsBins, GetTotalScalevars(), GetTotalScalenodes(), XmaxFromI[0], NSubproc);
         ResizeTable( &SigmaTilde , fNObsBins, GetTotalScalevars(), GetTotalScalenodes(), XmaxFromI, NSubproc );
         int nst = ReadTable  ( &SigmaTilde , table );
         //printf("  *  fastNLOCoefficients::Read(). Read %d lines of SigmaTilde.\n",nst);
         printf("  *  fastNLOCoefficients::Read(). Read %d lines of FASTNLO v2.0 tables.\n",nst+nsn);

      }


      if ( NScaleDep >= 3 ) {

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

         nn3 += ReadFlexibleVector  ( &ScaleNode1 , table );
         nn3 += ReadFlexibleVector  ( &ScaleNode2 , table );
         NscalenodeScale1 = ScaleNode1[0].size();
         NscalenodeScale2 = ScaleNode2[0].size();

         nn3 += ReadFlexibleVector  ( &SigmaTildeMuIndep , table , true );
         //if ( NScaleDep==3 || fScen->ILOord!=Npow || NScaleDep==5 ){
         if ( NScaleDep==3 || NScaleDep==5 ){
            //cout<<"Read NLO FlexTable. NScaleDep="<<NScaleDep<<"\tNpow="<<Npow<<"\tfScen->ILOord="<<fScen->ILOord<<endl;
            nn3 += ReadFlexibleVector  ( &SigmaTildeMuFDep , table , true );
            nn3 += ReadFlexibleVector  ( &SigmaTildeMuRDep , table , true );
         }
         nn3 += ReadFlexibleVector  ( &SigmaRefMixed , table , true );
         nn3 += ReadFlexibleVector  ( &SigmaRef_s1 , table , true );
         nn3 += ReadFlexibleVector  ( &SigmaRef_s2 , table , true );

//       *table >> NscalenodeScale1 ;
//       ResizeTable( &ScaleNode1 , fNObsBins , NscalenodeScale1 );
//       nn3 += ReadTable  ( &ScaleNode1 , table );

//       *table >> NscalenodeScale2 ;
//       ResizeTable( &ScaleNode2 , fNObsBins , NscalenodeScale2 );
//       nn3 += ReadTable  ( &ScaleNode2 , table );

//       int XMaxFromFromDim[1] = { 0 };
//       //ResizeTable( &PdfLcMuVar , fNObsBins , XMaxFromFromDim , NscalenodeScale1 , NscalenodeScale2 , NSubproc );
//       //ResizeTable( &AlphasTwoPi , fNObsBins , NscalenodeScale1 , NscalenodeScale2 );

//       ResizeTable( &SigmaTildeMuIndep , fNObsBins , XMaxFromFromDim , NscalenodeScale1 , NscalenodeScale2 , NSubproc );
//       nn3 += ReadTable  ( &SigmaTildeMuIndep , table );

//       ResizeTable( &SigmaTildeMuFDep , fNObsBins , XMaxFromFromDim , NscalenodeScale1 , NscalenodeScale2 , NSubproc );
//       nn3 += ReadTable  ( &SigmaTildeMuFDep , table );

//       ResizeTable( &SigmaTildeMuRDep , fNObsBins , XMaxFromFromDim , NscalenodeScale1 , NscalenodeScale2 , NSubproc );
//       nn3 += ReadTable  ( &SigmaTildeMuRDep , table );

//       ResizeTable( &SigmaRefMixed , fNObsBins , NSubproc );
//       nn3 += ReadTable  ( &SigmaRefMixed , table );

//       ResizeTable( &SigmaRef_s1 , fNObsBins , NSubproc );
//       nn3 += ReadTable  ( &SigmaRef_s1 , table );

//       ResizeTable( &SigmaRef_s2 , fNObsBins , NSubproc );
//       nn3 += ReadTable  ( &SigmaRef_s2 , table );
         printf("  *  fastNLOCoefficients::Read(). Read %d lines of NScaleDep>=3 Tables.\n",nn3);

      }


   }// end of not data and not corrections

   key = 0;
   *table >> key;
   if(key != tablemagicno){
      printf("fastNLOCoefficients::Read: At end of block found %d instead of %d.\n",key,tablemagicno);
      printf("                  You might have 'nan' in your table.\n");
      return 1;
   };
   // Put magic number back
   for(int i=0;i<(int)(log10((double)key)+1);i++){
      table->unget();
   }
   return 0;
}

int fastNLOCoefficients::Write(ostream *table, int option){
   *table << tablemagicno << sep;
   *table << IXsectUnits << sep;
   *table << IDataFlag << sep;
   *table << IAddMultFlag << sep;
   *table << IContrFlag1 << sep;
   *table << IContrFlag2 << sep;
   //KR: IContrFlag3 replaced by NScaleDep
   //*table << IContrFlag3 << sep;     // v2.0+. for v2.1 write IContrFlag3 here, but NScaleDep only later
   if ( NScaleDep==3 ) {
      if ( Npow==fILOord) {
         cout<<" * Increase NScaleDep from 3 to 4, because LO!"<<endl;
         NScaleDep=4;
      }
      else {
         cout<<" * Increase NScaleDep from 3 to 5 because NLO!"<<endl;
         NScaleDep=5;
      }
   }
   *table << NScaleDep << sep;
   *table << CtrbDescript.size() << sep;
   //printf("  *  fastNLOCoefficients::Write().  IDataFlag: %d, IAddMultFlag: %d, IContrFlag1: %d, IContrFlag2: %d, NScaleDep: %d\n",
   //IDataFlag,IAddMultFlag,IContrFlag1,IContrFlag2,NScaleDep);
   for(unsigned int i=0;i<CtrbDescript.size();i++){
      *table << CtrbDescript[i] << sep;
   }
   *table << CodeDescript.size() << sep;
   for(unsigned int i=0;i<CodeDescript.size();i++){
      *table << CodeDescript[i] << sep;
   }

   if(IDataFlag==1){
      *table << Nuncorrel << sep;
      for(int i=0;i<Nuncorrel;i++){
         *table << UncDescr[i] << sep;
      }
      *table << Ncorrel << sep;
      for(int i=0;i<Ncorrel;i++){
         *table << CorDescr[i]  << sep;
      }
      for(int i=0;i<fNObsBins;i++){
         *table << Xcenter[i] << sep;
         *table << Value[i] << sep;
         for(int j=0;j<Nuncorrel;j++){
            *table << UncorLo[i][j] << sep;
            *table << UncorHi[i][j] << sep;
         }
         for(int j=0;j<Ncorrel;j++){
            *table << CorrLo[i][j] << sep;
            *table << CorrHi[i][j] << sep;
         }
      }
      *table << NErrMatrix << sep;
      for(int i=0;i<NErrMatrix;i++){
         for(int j=0;j<(int)pow((double)fNObsBins,2);j++){
            *table << matrixelement[i][j] << sep;
         }
      }
   }// end of IDataFlag==1

   cout<<" 3"<<endl;
   if(IAddMultFlag==1){
      *table << Nuncorrel << sep;
      for(int i=0;i<Nuncorrel;i++){
         *table << UncDescr[i]  << sep;
      }
      *table << Ncorrel << sep;
      for(int i=0;i<Ncorrel;i++){
         *table << CorDescr[i]  << sep;
      }
      for(int i=0;i<fNObsBins;i++){
         *table << fact[i] << sep;
         for(int j=0;j<Nuncorrel;j++){
            *table << UncorLo[i][j] << sep;
            *table << UncorHi[i][j] << sep;
         }
         for(int j=0;j<Ncorrel;j++){
            *table << CorrLo[i][j] << sep;
            *table << CorrHi[i][j] << sep;
         }
      }
   }// end of IAddMultFlag==1

   if(!(IDataFlag==1) && !(IAddMultFlag==1)){
      *table << IRef << sep;
      *table << IScaleDep << sep;
      *table << Nevt << sep;
      *table << Npow << sep;
      *table << NPDF << sep;
      if(NPDF>0){
         for(int i=0;i<NPDF;i++){
            *table <<  NPDFPDG[i] << sep;
         }
      }
      *table << NPDFDim << sep;
      *table << NFragFunc << sep;
    if(NFragFunc>0){
         for(int i=0;i<NFragFunc;i++){
            *table <<  NFFPDG[i] << sep;
         }
      }
      *table << NFFDim << sep;
      *table << NSubproc << sep;
      *table << IPDFdef1 << sep;
      *table << IPDFdef2 << sep;
      *table << IPDFdef3 << sep;
     if(IPDFdef1==0){
         for(int i=0;i<NSubproc;i++){
            // Missing: linear PDF combinations for IPDFdef1=0
            if(NPDF==1){
            }else{
               if(NPDF==2){
               }
            }
         }
      }
      for(int i=0;i<fNObsBins;i++){
         *table << XNode1[i].size() << sep;
         for(unsigned int j=0;j<XNode1[i].size();j++){
            *table << XNode1[i][j] << sep;
         }
      }
      if(NPDFDim==2){
         for(int i=0;i<fNObsBins;i++){
            *table << XNode2[i].size() << sep;
            for(unsigned int j=0;j<XNode2[i].size();j++){
               *table << XNode2[i][j] << sep;
            }
         }
      }
      cout<<" 10"<<endl;
     if(NFragFunc>0){
         for(int i=0;i<fNObsBins;i++){
            *table << Nztot[i] << sep;
            for(int j=0;j<Nztot[i];j++){
               *table << ZNode[i][j] << sep;
            }
         }
      }
      *table << NScales << sep;
      *table << NScaleDim << sep;
      for(int i=0;i<NScales;i++){
         *table << Iscale[i] << sep;
      }
     for(int i=0;i<NScaleDim;i++){
         *table << ScaleDescript[i].size() << sep;
         for(unsigned int j=0;j<ScaleDescript[i].size();j++){
            *table << ScaleDescript[i][j] << sep;
         }
      }

      //! v2.1 store NScaleDep here
      //! *table << NScaleDep << sep;
      cout<<"fastNLOCoefficients. Writing coefficients."<<endl;

      if ( NScaleDep<3 ){
         for(int i=0;i<NScaleDim;i++){
            *table << Nscalevar[i] << sep;
            *table << Nscalenode[i] << sep;
         }
         for(int i=0;i<NScaleDim;i++){
            for(int j=0;j<Nscalevar[i];j++){
               *table << ScaleFac[i][j] << sep;
            }
         }

        int nsn = WriteTable( &ScaleNode  , table );
        //printf("  *  fastNLOCoefficients::Write(). Wrote %d lines of ScaleNode.\n",nsn);
        int nst = WriteTable( &SigmaTilde , table , (bool)(option & DividebyNevt) , Nevt );
        //printf("  *  fastNLOCoefficients::Write(). Wrote %d lines of SigmaTilde.\n",nst);
        printf("  *  fastNLOCoefficients::Write(). Wrote %d lines of FASTNLO v2.0 tables.\n",nst+nsn);


      } // end if NScaleDep !=3.

      if ( NScaleDep>=3 ) {
         int nn3 = 0;

         nn3 += WriteFlexibleTable( &ScaleNode1 , table );
         nn3 += WriteFlexibleTable( &ScaleNode2 , table );
         nn3 += WriteFlexibleTable( &SigmaTildeMuIndep, table , (bool)(option & DividebyNevt) , Nevt , true );
         //if ( NScaleDep==3 || Npow!=fScen->ILOord || NScaleDep==5) {
         if ( NScaleDep==3 || NScaleDep==5) {
            //cout<<"Write NLO FlexTable. NScaleDep="<<NScaleDep<<"\tNpow="<<Npow<<"\tfScen->ILOord="<<fScen->ILOord<<endl;
            nn3 += WriteFlexibleTable( &SigmaTildeMuFDep , table , (bool)(option & DividebyNevt) , Nevt , true );
            nn3 += WriteFlexibleTable( &SigmaTildeMuRDep , table , (bool)(option & DividebyNevt) , Nevt , true );
         }
         nn3 += WriteFlexibleTable( &SigmaRefMixed      , table , (bool)(option & DividebyNevt) , Nevt , true );
         nn3 += WriteFlexibleTable( &SigmaRef_s1        , table , (bool)(option & DividebyNevt) , Nevt , true );
         nn3 += WriteFlexibleTable( &SigmaRef_s2        , table , (bool)(option & DividebyNevt) , Nevt , true );

//       *table << NscalenodeScale1 << sep;
//       nn3 += WriteTable( &ScaleNode1 , table );

//       *table << NscalenodeScale2 << sep;
//       nn3 += WriteTable( &ScaleNode2 , table );

//       nn3 += WriteTable( &SigmaTildeMuIndep, table , (bool)(option & DividebyNevt) , Nevt );
//       nn3 += WriteTable( &SigmaTildeMuFDep , table , (bool)(option & DividebyNevt) , Nevt );
//       nn3 += WriteTable( &SigmaTildeMuRDep , table , (bool)(option & DividebyNevt) , Nevt );

//       nn3 += WriteTable( &SigmaRefMixed      , table , (bool)(option & DividebyNevt) , Nevt );
//       nn3 += WriteTable( &SigmaRef_s1        , table , (bool)(option & DividebyNevt) , Nevt );
//       nn3 += WriteTable( &SigmaRef_s2        , table , (bool)(option & DividebyNevt) , Nevt );

         printf("  *  fastNLOCoefficients::Write(). Wrote %d lines of v2.1 Tables.\n",nn3);

      } // if(NScaleDep==3)
   }// end of not data and not corrections

   return 0;
}

int fastNLOCoefficients::Copy(fastNLOCoefficients* other){

   streambuf* streambuf = new stringbuf(ios_base::in | ios_base::out);
   iostream* buffer = new iostream(streambuf);
   other->Write(buffer);
   *buffer << tablemagicno << endl;
   this->Read(buffer);
   delete buffer;
   delete streambuf;

   return(0);
}

void fastNLOCoefficients::Add(fastNLOCoefficients* other){
   double w1 = (double)Nevt / (Nevt+other->Nevt);
   double w2 = (double)other->Nevt / (Nevt+other->Nevt);
   Nevt += other->Nevt;

   if ( NScaleDep<3 ){
      AddTableToAnotherTable( &SigmaTilde , &(other->SigmaTilde) ,w1 , w2 );
   }

   if ( NScaleDep >= 3 ){
     AddTableToAnotherTable( &SigmaTildeMuIndep , &(other->SigmaTildeMuIndep) ,w1 , w2 );
     //if ( NScaleDep==3 || NScaleDep==5 || fScen->ILOord!=Npow) {
     if ( NScaleDep==3 || NScaleDep==5 ) {
        //cout<<"Adding NLO table"<<endl;
        AddTableToAnotherTable( &SigmaTildeMuFDep , &(other->SigmaTildeMuFDep) ,w1 , w2 );
        AddTableToAnotherTable( &SigmaTildeMuRDep , &(other->SigmaTildeMuRDep) ,w1 , w2 );
     }
     AddTableToAnotherTable( &SigmaRefMixed , &(other->SigmaRefMixed) ,w1 , w2 );
     AddTableToAnotherTable( &SigmaRef_s1 , &(other->SigmaRef_s1) ,w1 , w2 );
     AddTableToAnotherTable( &SigmaRef_s2 , &(other->SigmaRef_s2) ,w1 , w2 );
   }

}

void fastNLOCoefficients::StripWhitespace(string &str) const {
   for(string::iterator achar = str.end(); achar>str.begin();achar--) {
      if (*achar==0x20 || *achar==0x00){
         str.erase(achar);
      }else{
         break;
      }
   }
}

int fastNLOCoefficients::GetNxmax(int i) const {
   int nxmax = 0;
   switch (NPDFDim) {
   case 0: nxmax = (int)XNode1[i].size();
      break;
      //   case 1: nxmax = ((int)pow((double)Nxtot1[i],2)+Nxtot1[i])/2;
   case 1: nxmax = ((int)pow((double)XNode1[i].size(),2)+XNode1[i].size())/2;
      break;
   case 2: nxmax = XNode1[i].size()*XNode2[i].size();
      break;
   default: ;
   }
   return nxmax;
};

int fastNLOCoefficients::GetXIndex(int Obsbin,int x1bin,int x2bin) const {
   int xbin = 0;
   switch (NPDFDim) {
   case 0: xbin = x1bin; // linear
      break;
   case 1: xbin = x1bin + (x2bin*(x2bin+1)/2);    // half matrix
      break;
   case 2: xbin = x1bin + x2bin * XNode1[Obsbin].size(); // full matrix
      break;
   default: ;
   }
   return xbin;
};



int fastNLOCoefficients::GetTotalScalevars() const {
   int totalscalevars=1;
   for(int scaledim=0;scaledim<NScaleDim;scaledim++){
      totalscalevars *= Nscalevar[scaledim];
   }
   return totalscalevars;
}

int fastNLOCoefficients::GetTotalScalenodes() const {
   int totalscalenodes=1;
   for(int scaledim=0;scaledim<NScaleDim;scaledim++){
      totalscalenodes *= Nscalenode[scaledim];
   }
   return totalscalenodes;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoefficients::ResizeTable(vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int* dim5GetNxmaxFromDimI , int dim6 ){
  if ( dim0 > 0 ){
    if ( dim5GetNxmaxFromDimI[0] == 0 ) {
      v->resize(dim0);
      for ( int i= 0 ; i<dim0 ; i++){
        ResizeTable( &(v->at(i)) , dim1, dim2, dim3, dim4, GetNxmax(i), dim6 );
      }
    }
    else if ( dim5GetNxmaxFromDimI[0] != 0 ){
      cout << "Error in Resize Table. This is not yet implemented" << endl;
      exit(1);
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }

}


//________________________________________________________________________________________________________________ //
void fastNLOCoefficients::ResizeTable(vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5, int dim6 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
    for ( int i= 0 ; i<dim0 ; i++){
      ResizeTable( &(v->at(i)) , dim1, dim2, dim3, dim4, dim5, dim6 );
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoefficients::ResizeTable( vector<vector<vector<vector<vector<vector<double > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
    for ( int i= 0 ; i<dim0 ; i++){
      ResizeTable( &(v->at(i)) , dim1, dim2, dim3, dim4, dim5 );
    }
  } else {
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoefficients::ResizeTable( vector<vector<vector<vector<vector<double > > > > >* v, int dim0 , int dim1, int dim2, int* dim3GetNxmaxFromDimI, int dim4 ){
  if ( dim0 > 0 ){
    if ( dim3GetNxmaxFromDimI[0] == 0 ) {
      v->resize(dim0);
      for ( int i= 0 ; i<dim0 ; i++){
        ResizeTable( &(v->at(i)) , dim1, dim2, GetNxmax(i), dim4 );
      }
    }
    else if ( dim3GetNxmaxFromDimI[0] != 0 ){
      cout << "Error in Resize Table. This is not yet implemented" << endl;
      exit(1);
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoefficients::ResizeTable( vector<vector<vector<vector<vector<double > > > > >* v, int dim0 , int* dim1GetNxmaxFromDimI, int dim2, int dim3, int dim4 ){
  if ( dim0 > 0 ){
    if ( dim1GetNxmaxFromDimI[0] == 0 ) {
      v->resize(dim0);
      for ( int i= 0 ; i<dim0 ; i++){
        ResizeTable( &(v->at(i)) , GetNxmax(i), dim2, dim3, dim4 );
      }
    }
    else if ( dim1GetNxmaxFromDimI[0] != 0 ){
      cout << "Error in Resize Table. This is not yet implemented" << endl;
      exit(1);
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }

}


//________________________________________________________________________________________________________________ //
void fastNLOCoefficients::ResizeTable( vector<vector<vector<vector<vector<double > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
    for ( int i= 0 ; i<dim0 ; i++){
      ResizeTable( &(v->at(i)) , dim1, dim2, dim3, dim4 );
    }
  } else {
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoefficients::ResizeTable( vector<vector<vector<vector<double > > > >* v, int dim0 , int dim1, int* dim2GetNxmaxFromDimI, int dim3 ){
  if ( dim0 > 0 ){
    if ( dim2GetNxmaxFromDimI[0] == 0 ) {
      v->resize(dim0);
      for ( int i= 0 ; i<dim0 ; i++){
        ResizeTable( &(v->at(i)) , dim1, GetNxmax(i), dim3 );
      }
    }
    else if ( dim2GetNxmaxFromDimI[0] != 0 ){
      cout << "Error in Resize Table. This is not yet implemented" << endl;
      exit(1);
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }

}


//________________________________________________________________________________________________________________ //
void fastNLOCoefficients::ResizeTable( vector<vector<vector<vector<double > > > >* v, int dim0 , int dim1, int dim2, int dim3 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
    for ( int i= 0 ; i<dim0 ; i++){
      ResizeTable( &(v->at(i)) , dim1, dim2, dim3 );
    }
  } else {
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoefficients::ResizeTable( vector<vector<vector<double > > >* v, int dim0 , int* dim1GetNxmaxFromDimI, int dim2 ){
  if ( dim0 > 0 ){
    if ( dim1GetNxmaxFromDimI[0] == 0 ) {
      v->resize(dim0);
      for ( int i= 0 ; i<dim0 ; i++){
        ResizeTable( &(v->at(i)) , GetNxmax(i), dim2 );
      }
    }
    else if ( dim1GetNxmaxFromDimI[0] != 0 ){
      cout << "Error in Resize Table. This is not yet implemented" << endl;
      exit(1);
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }

}


//________________________________________________________________________________________________________________ //
void fastNLOCoefficients::ResizeTable( vector<vector<vector<double > > >* v, int dim0 , int dim1, int dim2 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
    for ( int i= 0 ; i<dim0 ; i++){
      ResizeTable( &(v->at(i)) , dim1, dim2 );
    }
  } else {
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoefficients::ResizeTable( vector<vector<double > >*  v, int dim0 , int* dim1GetNxmaxFromDimI ){
  if ( dim0 > 0 ){
    if ( dim1GetNxmaxFromDimI[0] == 0 ) {
      v->resize(dim0);
      for ( int i= 0 ; i<dim0 ; i++){
        ResizeTable( &(v->at(i)) , GetNxmax(i) );
      }
    }
    else if ( dim1GetNxmaxFromDimI[0] != 0 ){
      cout << "Error in Resize Table. This is not yet implemented" << endl;
      exit(1);
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoefficients::ResizeTable( vector<vector<double > >*  v, int dim0 , int dim1 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
    for ( int i= 0 ; i<dim0 ; i++){
      ResizeTable( &(v->at(i)) , dim1 );
    }
  } else {
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoefficients::ResizeTable( vector<double >* v, int dim0 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
  }
  else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}



//________________________________________________________________________________________________________________ //
int fastNLOCoefficients::ReadTable(vector<double >* v, istream *table ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    *table >> v->at(i0);
    nn++;
  }
  return nn;
}


//________________________________________________________________________________________________________________ //
int fastNLOCoefficients::WriteTable(vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, ostream *table , bool DivByNevt , int Nevt  ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
        for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
          for(unsigned int i4=0;i4<v->at(i0)[i1][i2][i3].size();i4++){
            for(unsigned int i5=0;i5<v->at(i0)[i1][i2][i3][i4].size();i5++){
              for(unsigned int i6=0;i6<v->at(i0)[i1][i2][i3][i4][i5].size();i6++){
                if( DivByNevt && Nevt>0){
                  *table << v->at(i0)[i1][i2][i3][i4][i5][i6] / Nevt << sep;
                }else{
                  *table << v->at(i0)[i1][i2][i3][i4][i5][i6] << sep;
                }
                nn++;
              }
            }
          }
        }
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fastNLOCoefficients::WriteTable(vector<vector<vector<vector<vector<vector<double > > > > > >* v, ostream *table , bool DivByNevt , int Nevt ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
        for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
          for(unsigned int i4=0;i4<v->at(i0)[i1][i2][i3].size();i4++){
            for(unsigned int i5=0;i5<v->at(i0)[i1][i2][i3][i4].size();i5++){
              if( DivByNevt && Nevt>0){
                *table << v->at(i0)[i1][i2][i3][i4][i5] / Nevt << sep;
              }else{
                *table << v->at(i0)[i1][i2][i3][i4][i5] << sep;
              }
              nn++;
            }
          }
        }
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fastNLOCoefficients::WriteTable(vector<vector<vector<vector<vector<double > > > > >* v, ostream *table , bool DivByNevt , int Nevt ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
        for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
          for(unsigned int i4=0;i4<v->at(i0)[i1][i2][i3].size();i4++){
            if( DivByNevt && Nevt>0){
              *table << v->at(i0)[i1][i2][i3][i4] / Nevt << sep;
                }else{
              *table << v->at(i0)[i1][i2][i3][i4] << sep;
            }
            nn++;
          }
        }
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fastNLOCoefficients::WriteTable(vector<vector<vector<vector<double > > > >* v, ostream *table , bool DivByNevt , int Nevt ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
        for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
          if( DivByNevt && Nevt>0){
            *table << v->at(i0)[i1][i2][i3] / Nevt << sep;
          }else{
            *table << v->at(i0)[i1][i2][i3] << sep;
          }
          nn++;
        }
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fastNLOCoefficients::WriteTable(vector<vector<vector<double > > >* v, ostream *table , bool DivByNevt , int Nevt ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
        if( DivByNevt && Nevt>0){
          *table << v->at(i0)[i1][i2] / Nevt << sep;
        }else{
          *table << v->at(i0)[i1][i2] << sep;
        }
        nn++;
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fastNLOCoefficients::WriteTable(vector<vector<double > >* v, ostream *table , bool DivByNevt , int Nevt ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      if( DivByNevt && Nevt>0){
        *table << v->at(i0)[i1] / Nevt << sep;
      }else{
        *table << v->at(i0)[i1] << sep;
      }
      nn++;
    }
  }
  return nn;
}



//________________________________________________________________________________________________________________ //
int fastNLOCoefficients::WriteTable(vector<double >* v, ostream *table , bool DivByNevt , int Nevt ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    if( DivByNevt && Nevt>0){
      *table << v->at(i0) / Nevt << sep;
    }else{
      *table << v->at(i0) << sep;
    }
    nn++;
  }
  return nn;
}


//________________________________________________________________________________________________________________ //
int fastNLOCoefficients::WriteFlexibleTable(vector<vector<vector<vector<vector<vector<vector< double > > > > > > >* v, ostream *table , bool DivByNevt , int Nevt , bool nProcLast ){
   int nn = 1;
   *table << v->size() << sep;
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += WriteFlexibleTable( &(v->at(i0)) , table , DivByNevt, Nevt , nProcLast );
   }
   return nn;
}
//________________________________________________________________________________________________________________ //
int fastNLOCoefficients::WriteFlexibleTable(vector<vector<vector<vector<vector<vector<double > > > > > >* v, ostream *table , bool DivByNevt , int Nevt , bool nProcLast ){
   int nn = 1;
   *table << v->size() << sep;
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += WriteFlexibleTable( &(v->at(i0)) , table , DivByNevt, Nevt , nProcLast );
   }
   return nn;
}
//________________________________________________________________________________________________________________ //
int fastNLOCoefficients::WriteFlexibleTable(vector<vector<vector<vector<vector<double > > > >  >* v, ostream *table , bool DivByNevt , int Nevt , bool nProcLast ){
   int nn = 1;
   *table << v->size() << sep;
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += WriteFlexibleTable( &(v->at(i0)) , table , DivByNevt, Nevt , nProcLast );
   }
   return nn;
}
//________________________________________________________________________________________________________________ //
int fastNLOCoefficients::WriteFlexibleTable(vector<vector<vector<vector<double > >  > >* v, ostream *table , bool DivByNevt , int Nevt , bool nProcLast ){
   int nn = 1;
   *table << v->size() << sep;
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += WriteFlexibleTable( &(v->at(i0)) , table , DivByNevt, Nevt , nProcLast );
   }
   return nn;
}
//________________________________________________________________________________________________________________ //
int fastNLOCoefficients::WriteFlexibleTable(vector<vector<vector<double > > >* v, ostream *table , bool DivByNevt , int Nevt , bool nProcLast ){
   int nn = 1;
   *table << v->size() << sep;
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += WriteFlexibleTable( &(v->at(i0)) , table , DivByNevt, Nevt , nProcLast );
   }
   return nn;
}
//________________________________________________________________________________________________________________ //
int fastNLOCoefficients::WriteFlexibleTable(vector<vector<double > >* v, ostream *table , bool DivByNevt , int Nevt , bool nProcLast ){
   int nn = 1;
   *table << v->size() << sep;
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += WriteFlexibleTable( &(v->at(i0)) , table , DivByNevt, Nevt , nProcLast );
   }
   return nn;
}

//________________________________________________________________________________________________________________ //
int fastNLOCoefficients::WriteFlexibleTable(vector<double >* v, ostream *table , bool DivByNevt , int Nevt , bool nProcLast ){
   int nn = 1;
   if ( !nProcLast )*table << v->size() << sep;
   for(unsigned int i0=0;i0<v->size();i0++){
      if( DivByNevt && Nevt>0)  *table << v->at(i0) / Nevt << sep;
      else                      *table << v->at(i0) << sep;
      nn++;
   }
   return nn;
}

//________________________________________________________________________________________________________________ //
int fastNLOCoefficients::ReadFlexibleVector(vector<double >* v, istream *table , bool nProcLast ){
   int nn = 0;
   if ( !nProcLast ) {
      int size = 0;
      *table >> size; nn++;
      v->resize(size);
   }
   else {
      v->resize(NSubproc);
   }
   for(unsigned int i0=0;i0<v->size();i0++){
      *table >> v->at(i0);
      nn++;
   }
   return nn;
}



//________________________________________________________________________________________________________________ //
void fastNLOCoefficients::AddTableToAnotherTable( vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* vSum, vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* vAdd, double w1, double w2){
  if ( vSum->size() != vAdd->size() ) {cout<<"Error in fastNLOCoefficients::AddTableToAnotherTable. Cannot add tables with different size. [v7] s1="<<vSum->size()<<", s2="<<vAdd->size()<<endl; return;}
  for ( unsigned int i = 0 ; i<vSum->size() ; i++ ){
    AddTableToAnotherTable( &(vSum->at(i)), &(vAdd->at(i)), w1 , w2  );
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoefficients::AddTableToAnotherTable( vector<vector<vector<vector<vector<vector<double > > > > > >* vSum, vector<vector<vector<vector<vector<vector<double > > > > > >* vAdd, double w1, double w2){

  if ( vSum->size() != vAdd->size() ) {cout<<"Error in fastNLOCoefficients::AddTableToAnotherTable. Cannot add tables with different size. [v6] s1="<<vSum->size()<<", s2="<<vAdd->size()<<endl; return;}
  for ( unsigned int i = 0 ; i<vSum->size() ; i++ ){
    AddTableToAnotherTable( &(vSum->at(i)), &(vAdd->at(i)), w1 , w2  );
  }
}

//________________________________________________________________________________________________________________ //


void fastNLOCoefficients::AddTableToAnotherTable( vector<vector<vector<vector<vector<double > > > > >* vSum, vector<vector<vector<vector<vector<double > > > > >* vAdd, double w1, double w2){
  if ( vSum->size() != vAdd->size() ) {cout<<"Error in fastNLOCoefficients::AddTableToAnotherTable. Cannot add tables with different size. [v5] s1="<<vSum->size()<<", s2="<<vAdd->size()<<endl; return;}
  for ( unsigned int i = 0 ; i<vSum->size() ; i++ ){
    AddTableToAnotherTable( &(vSum->at(i)), &(vAdd->at(i)), w1 , w2  );
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoefficients::AddTableToAnotherTable( vector<vector<vector<vector<double > > > >* vSum, vector<vector<vector<vector<double > > > >* vAdd, double w1, double w2){
  if ( vSum->size() != vAdd->size() ) {cout<<"Error in fastNLOCoefficients::AddTableToAnotherTable. Cannot add tables with different size. [v4] s1="<<vSum->size()<<", s2="<<vAdd->size()<<endl; return;}
  for ( unsigned int i = 0 ; i<vSum->size() ; i++ ){
    AddTableToAnotherTable( &(vSum->at(i)), &(vAdd->at(i)), w1 , w2  );
  }
}

//________________________________________________________________________________________________________________ //

void fastNLOCoefficients::AddTableToAnotherTable( vector<vector<vector<double > > >* vSum, vector<vector<vector<double > > >* vAdd, double w1, double w2){

  if ( vSum->size() != vAdd->size() ) {cout<<"Error in fastNLOCoefficients::AddTableToAnotherTable. Cannot add tables with different size. [v3] s1="<<vSum->size()<<", s2="<<vAdd->size()<<endl; return;}
  for ( unsigned int i = 0 ; i<vSum->size() ; i++ ){
    AddTableToAnotherTable( &(vSum->at(i)), &(vAdd->at(i)), w1 , w2  );
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoefficients::AddTableToAnotherTable( vector<vector<double > >* vSum, vector<vector<double > >* vAdd, double w1, double w2){
  if ( vSum->size() != vAdd->size() ) {cout<<"Error in fastNLOCoefficients::AddTableToAnotherTable. Cannot add tables with different size. [v2] s1="<<vSum->size()<<", s2="<<vAdd->size()<<endl; return;}
  for ( unsigned int i = 0 ; i<vSum->size() ; i++ ){
    AddTableToAnotherTable( &(vSum->at(i)), &(vAdd->at(i)), w1 , w2  );
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoefficients::AddTableToAnotherTable( vector<double >* vSum, vector<double >* vAdd, double w1, double w2){
   if ( vSum->size() != vAdd->size() ) {cout<<"Error in fastNLOCoefficients::AddTableToAnotherTable. Cannot add tables with different size. [v1] s1="<<vSum->size()<<", s2="<<vAdd->size()<<endl; return;}
  for ( unsigned int i = 0 ; i<vSum->size() ; i++ ){
    (*vSum)[i] =  w1*(*vSum)[i] + w2*(*vAdd)[i];
  }
}



//________________________________________________________________________________________________________________ //


void fastNLOCoefficients::SetNlojetDefaults(){
  SetIDataFlag(0);
  SetIAddMultFlag(0);
  SetIContrFlag1(1);
  SetIContrFlag2(100); // specify if LO or NLO
  SetNScaleDep(0);
  SetIXsectUnits(12);
  SetNlojetDescr();
};


void fastNLOCoefficients::Print() const {
  printf("\n **************** FastNLO Table: BlockB ****************\n\n");
  printf(" B   Scenario::GetNObsBin()        %d\n",fNObsBins);
  printf(" B   IXsectUnits                   %d\n",IXsectUnits);
  printf(" B   IDataFlag                     %d\n",IDataFlag);
  printf(" B   IAddMultFlag                  %d\n",IAddMultFlag);
  printf(" B   IContrFlag1                   %d\n",IContrFlag1);
  printf(" B   IContrFlag2                   %d\n",IContrFlag2);
  //  printf(" B   IContrFlag3 (always 0)        %d\n",IContrFlag3);
  printf(" B   NScaleDep                     %d\n",NScaleDep);
  for(unsigned int i=0;i<CtrbDescript.size();i++){
    printf(" B   CtrbDescript[%d]               %s\n",i,CtrbDescript[i].data());
  }
  //printf(" B   NCodeDescr                    %d\n",NCodeDescr);
  for(unsigned int i=0;i<CodeDescript.size();i++){
    printf(" B   CodeDescript[%d]               %s\n",i,CodeDescript[i].data());
  }

  if(IDataFlag==1){
    printf(" B   Nuncorrel                     %d\n",Nuncorrel);
    printf(" B   Ncorrel                       %d\n",Ncorrel);
    printf(" B   NErrMatrix                    %d\n",NErrMatrix);
    printf(" B   some more output could be printed here (IDataFlag==1).\n");
  }
  if(IAddMultFlag==1){
    printf(" B   some more output could be printed here (IAddMultFlag==1).\n");
  }

  if(!(IDataFlag==1) && !(IAddMultFlag==1)){ // that's the usual case
    printf(" B   IRef                          %d\n",IRef);
    printf(" B   IScaleDep                     %d\n",IScaleDep);
    printf(" B   Nevt                          %llu\n",Nevt);
    printf(" B   Npow                          %d\n",Npow);
    printf(" B   NPDF                          %d\n",NPDF);
    if(NPDF>0){
      for(int i=0;i<NPDF;i++){
        printf(" B    - NPDFPDG[%d]                 %d\n",i,NPDFPDG[i]);
      }
    }
    printf(" B   NPDFDim                       %d\n",NPDFDim);
    printf(" B   NFragFunc                     %d\n",NFragFunc);
    if(NFragFunc>0){
      for(int i=0;i<NFragFunc;i++){
        printf(" B    - NFFPDG[%d]               %d\n",i,NFFPDG[i]);
      }
    }
    printf(" B   NFFDim                        %d\n",NFFDim);
    printf(" B   NSubproc                      %d\n",NSubproc);
    printf(" B   IPDFdef1                      %d\n",IPDFdef1);
    printf(" B   IPDFdef2                      %d\n",IPDFdef2);
    printf(" B   IPDFdef3                      %d\n",IPDFdef3);
    printf(" B   Nxtot1[0-%d]             ",fNObsBins);
    for(int i=0;i<fNObsBins;i++){
       printf("%lu ,",XNode1[i].size());
    }
    printf(" B   \n");

//     for(int i=0;i<fNObsBins;i++){
//       printf(" B    XNode1[%d]             ",i);
//       for(int j=0;j<Nxtot1[i];j++){
//      printf(" B   %8.4f ,",XNode1[i][j]);
//       }
//       printf(" B   \n");
//     }
    printf(" B   if (NPDFDim==2), you could print xnodes2 here. (NPDFDim = %d)\n",NPDFDim);
    printf(" B   if (NFragFunc>0), you could print xnodes2 here. (NFragFunc = %d)\n",NFragFunc);
    printf(" B   NScales                       %d\n",NScales);
    for(int i=0;i<NScales;i++){
      printf(" B    - Iscale[%d]                  %d\n",i,Iscale[i]);
    }
    printf(" B   NScaleDim                     %d\n",NScaleDim);
    for(int i=0;i<NScaleDim;i++){
       //printf(" B    -  NscaleDescript[%d]         %d\n",i,NscaleDescript[i]);
       for(unsigned int j=0;j<ScaleDescript[i].size();j++){
        printf(" B    -  - ScaleDescript[%d][%d]     %s\n",i,j,ScaleDescript[i][j].data());
      }
       if ( NScaleDep<3 ) {
          printf(" B    - Nscalenode[%d]              %d\n",i,Nscalenode[i]);
          printf(" B    - Nscalevar[%d]               %d\n",i,Nscalevar[i]);
          for(int j=0;j<Nscalevar[i];j++){
             printf(" B    -  - ScaleFac[%d][%d]          %6.4f\n",i,j,ScaleFac[i][j]);
          }
       }
    }
    printf(" B   No printing of ScaleNode implemented yet.\n");
    printf(" B   No printing of SigmaTilde implemented yet.\n");
    if ( NScaleDep == 2 )
      printf(" B   NScaleDep == 2 :              yes.\n");
    if ( NScaleDep == 2 ) {
      printf(" B   No printing of SigmaTilde2Scales, and Scale2Nodes, etc... implemented yet.\n");
    }
    if ( NScaleDep>=3 ) {
      printf(" B   NscalenodeScale1              %d\n",NscalenodeScale1);
      printf(" B   NscalenodeScale2              %d\n",NscalenodeScale2);
    }

  }
  printf("\n *******************************************************\n\n");


}


//________________________________________________________________________________________________________________ //
