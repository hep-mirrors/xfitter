// Author: Daniel Britzger
// DESY, 23/07/2011

//  Version 0.2, 

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  Data storage class for 'BlockB'-variables                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "FastNLOBlockB.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace std;


//______________________________________________________________________________



FastNLOBlockB::FastNLOBlockB(const char* name , const int NObsBins )
{
  fname		= (char*)name;
  fNObsBins	= NObsBins;
}


FastNLOBlockB::FastNLOBlockB(const char* name , const int NObsBins , istream* table)
{
  fname		= (char*)name;
  fNObsBins	= NObsBins;
  ReadBlockB(table);
}

FastNLOBlockB::~FastNLOBlockB(void)
{
}
//______________________________________________________________________________


void FastNLOBlockB::ReadBlockB(istream *table){

   table->peek();
   if (table->eof()){
      printf("FastNLOBlockB::Read: Cannot read from file.\n");
      return;
   }
   
   int key = 0;
   *table >> key;
   if(key != tablemagicno){
     printf("FastNLOBlockB::Read: At beginning of block found %d instead of %d.\n",key,tablemagicno);
     return;
   };

   *table >> IXsectUnits;
   *table >> IDataFlag;
   *table >> IAddMultFlag;
   *table >> IContrFlag1;
   *table >> IContrFlag2;
   IContrFlag3 = 0;
   *table >> NScaleDep;
   *table >> NContrDescr;
   //   printf("  *  FastNLOBlockB::Read().  IDataFlag: %d, IAddMultFlag: %d, IContrFlag1: %d, IContrFlag2: %d, IContrFlag3: %d, NScaleDep: %d\n",IDataFlag,IAddMultFlag,IContrFlag1,IContrFlag2,IContrFlag3,NScaleDep );
   CtrbDescript.resize(NContrDescr);
   char buffer[257];
   table->getline(buffer,256);
   for(int i=0;i<NContrDescr;i++){
      table->getline(buffer,256);
      CtrbDescript[i] = buffer;
      //      StripWhitespace(CtrbDescript[i]);
   }

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
      //printf("  *  FastNLOBlockB::Read(). IRef : %d, IScaleDep: %d, Nevt: %d, Npow: %d, NPDF: %d, NPDFDim: %d\n", IRef ,IScaleDep  ,Nevt  , Npow ,NPDF , NPDFDim  );

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
      Nxtot1.resize(fNObsBins);
      XNode1.resize(fNObsBins);
      for(int i=0;i<fNObsBins;i++){
         *table >> Nxtot1[i];
         XNode1[i].resize(Nxtot1[i]);
         for(int j=0;j<Nxtot1[i];j++){
            *table >> XNode1[i][j];
         }
      }
      if(NPDFDim==2){
         Nxtot2.resize(fNObsBins);
         XNode2.resize(fNObsBins);
         for(int i=0;i<fNObsBins;i++){
            *table >> Nxtot2[i];
            XNode2[i].resize(Nxtot2[i]);
            for(int j=0;j<Nxtot2[i];j++){
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
      NscaleDescript.resize(NScaleDim);
      ScaleDescript.resize(NScaleDim);
      for(int i=0;i<NScaleDim;i++){
         *table >> NscaleDescript[i];
         ScaleDescript[i].resize(NscaleDescript[i]);
         table->getline(buffer,256);
         for(int j=0;j<NscaleDescript[i];j++){
            table->getline(buffer,256);
            ScaleDescript[i][j] = buffer;
            //            StripWhitespace(ScaleDescript[i][j]);
         }
      }

      Nscalevar.resize(NScaleDim);
      Nscalenode.resize(NScaleDim);
      for(int i=0;i<NScaleDim;i++){
         *table >> Nscalevar[i];
         *table >> Nscalenode[i];
      }
      //printf("  *  FastNLOBlockB::Read().bins %d, NScalevar[0] %d, Nscalenode[0] %d,  NScaleDim %d  \n",fNObsBins, Nscalevar[0] , Nscalenode[0] , NScaleDim );

      ScaleFac.resize(NScaleDim);
      for(int i=0;i<NScaleDim;i++){
         ScaleFac[i].resize(Nscalevar[i]);
         for(int j=0;j<Nscalevar[i];j++){
            *table >> ScaleFac[i][j];
         }
      }

      ResizeTable( &ScaleNode , fNObsBins, 1 , Nscalevar[0] , Nscalenode[0] ); // should work, since NScaleDim==1 
      ReadTable  ( &ScaleNode , table );
      //printf("  *  FastNLOBlockB::Read(). Read %d lines of ScaleNode.\n",nsn);

      int XmaxFromI[1] = {0};
      ResizeTable( &SigmaTilde , fNObsBins, GetTotalScalevars(), GetTotalScalenodes(), XmaxFromI, NSubproc );
      ReadTable  ( &SigmaTilde , table );
      //printf("  *  FastNLOBlockB::Read(). Read %d lines of SigmaTilde.\n",nst);

      ResizeTable( &PdfLc , fNObsBins, GetTotalScalenodes(), XmaxFromI, NSubproc );
      ResizeTable( &AlphasTwoPi_v20 , fNObsBins, GetTotalScalenodes() );

      if ( NScaleDep == 3 ) {
        int nn3 = 0;

        *table >> NscalenodeScaleQ ;
        ResizeTable( &ScaleNodeQ , fNObsBins , NscalenodeScaleQ );
        nn3 += ReadTable  ( &ScaleNodeQ , table );

        *table >> NscalenodeScalePt ;
        ResizeTable( &ScaleNodePt , fNObsBins , NscalenodeScalePt );
        nn3 += ReadTable  ( &ScaleNodePt , table );

        int XMaxFromFromDim[1] = { 0 };
        ResizeTable( &PdfLcMuVar , fNObsBins , XMaxFromFromDim , NscalenodeScaleQ , NscalenodeScalePt , NSubproc );
        ResizeTable( &AlphasTwoPi , fNObsBins , NscalenodeScaleQ , NscalenodeScalePt );

        ResizeTable( &SigmaTildeMuIndep , fNObsBins , XMaxFromFromDim , NscalenodeScaleQ , NscalenodeScalePt , NSubproc );
        nn3 += ReadTable  ( &SigmaTildeMuIndep , table );

        ResizeTable( &SigmaTildeMuFDep , fNObsBins , XMaxFromFromDim , NscalenodeScaleQ , NscalenodeScalePt , NSubproc );
        nn3 += ReadTable  ( &SigmaTildeMuFDep , table );

        ResizeTable( &SigmaTildeMuRDep , fNObsBins , XMaxFromFromDim , NscalenodeScaleQ , NscalenodeScalePt , NSubproc );
        nn3 += ReadTable  ( &SigmaTildeMuRDep , table );

        ResizeTable( &SigmaRefMixed , fNObsBins , NSubproc );
        nn3 += ReadTable  ( &SigmaRefMixed , table );

        ResizeTable( &SigmaRefQ2 , fNObsBins , NSubproc );
        nn3 += ReadTable  ( &SigmaRefQ2 , table );

        ResizeTable( &SigmaRefMufQ2MuRMixed , fNObsBins , NSubproc );
        nn3 += ReadTable  ( &SigmaRefMufQ2MuRMixed , table );
        //printf(" *  FastNLOBlockB::Read(). Read %d lines of NScaleDep==3 Tables.\n",nn3);

      }

   }// end of not data and not corrections

   key = 0;
   *table >> key;
   if(key != tablemagicno){
      printf("FastNLOBlockB::Read: At end of block found %d instead of %d.\n",key,tablemagicno);
      return;
   };
   // Put magic number back
   for(int i=0;i<(int)(log10((double)key)+1);i++){
      table->unget();
   }

   //printf("... FastNLOBlockB::Reading: Table read succesfully .\n");

   return;
}


//______________________________________________________________________________


void FastNLOBlockB::Print(){
  printf("\n **************** BlockB: %s ****************\n\n",fname);
  printf(" B   fNNObsBins		     %d\n",fNObsBins);
  printf(" B   IXsectUnits                   %d\n",IXsectUnits);
  printf(" B   IDataFlag                     %d\n",IDataFlag);
  printf(" B   IAddMultFlag                  %d\n",IAddMultFlag);
  printf(" B   IContrFlag1                   %d\n",IContrFlag1);
  printf(" B   IContrFlag2                   %d\n",IContrFlag2);
  printf(" B   IContrFlag3 (always 0)        %d\n",IContrFlag3);
  printf(" B   NScaleDep                     %d\n",NScaleDep);
  printf(" B   NContrDescr                   %d\n",NContrDescr);
  for(int i=0;i<NContrDescr;i++){
    printf(" B   CtrbDescript[%d]               %s\n",i,CtrbDescript[i].data());
  }
  printf(" B   NCodeDescr                    %d\n",NCodeDescr);
  for(int i=0;i<NCodeDescr;i++){
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
    printf(" B   Nevt                          %e\n",Nevt*1.);
    printf(" B   Nevt                          %5.3e\n",Nevt*1.);
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
      printf("%d ,",Nxtot1[i]);
    }
    printf(" B   \n");
    printf(" B   if (NPDFDim==2), you could print xnodes2 here. (NPDFDim = %d)\n",NPDFDim);
    printf(" B   if (NFragFunc>0), you could print xnodes2 here. (NFragFunc = %d)\n",NFragFunc);
    printf(" B   NScales                       %d\n",NScales);
    for(int i=0;i<NScales;i++){
      printf(" B    - Iscale[%d]                  %d\n",i,Iscale[i]);
    }
    printf(" B   NScaleDim                     %d\n",NScaleDim);
    for(int i=0;i<NScaleDim;i++){
      printf(" B    -  NscaleDescript[%d]         %d\n",i,NscaleDescript[i]);
      for(int j=0;j<NscaleDescript[i];j++){
        printf(" B    -  - ScaleDescript[%d][%d]     %s\n",i,j,ScaleDescript[i][j].data());
      }
      printf(" B    - Nscalenode[%d]              %d\n",i,Nscalenode[i]);
      printf(" B    - Nscalevar[%d]               %d\n",i,Nscalevar[i]);
      for(int j=0;j<Nscalevar[i];j++){
        printf(" B    -  - ScaleFac[%d][%d]          %6.4f\n",i,j,ScaleFac[i][j]);
      }
    }
    printf(" B   No printing of ScaleNode implemented yet.\n");
    printf(" B   No printing of SigmaTilde implemented yet.\n");
    if ( NScaleDep == 3 ) {
      printf(" B   NscalenodeScaleQ              %d\n",NscalenodeScaleQ);
      printf(" B   NscalenodeScalePt             %d\n",NscalenodeScalePt);
    }

  }
  printf("\n *******************************************************\n\n");

}



//______________________________________________________________________________


//________________________________________________________________________________________________________________ //
int FastNLOBlockB::ReadTable(vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, istream *table ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
 	for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
	  for(unsigned int i4=0;i4<v->at(i0)[i1][i2][i3].size();i4++){
	    for(unsigned int i5=0;i5<v->at(i0)[i1][i2][i3][i4].size();i5++){
	      for(unsigned int i6=0;i6<v->at(i0)[i1][i2][i3][i4][i5].size();i6++){
		*table >> v->at(i0)[i1][i2][i3][i4][i5][i6];
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
int FastNLOBlockB::ReadTable(vector<vector<vector<vector<vector<vector<double > > > > > >* v, istream *table ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
 	for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
	  for(unsigned int i4=0;i4<v->at(i0)[i1][i2][i3].size();i4++){
	    for(unsigned int i5=0;i5<v->at(i0)[i1][i2][i3][i4].size();i5++){
	      *table >> v->at(i0)[i1][i2][i3][i4][i5];
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
int FastNLOBlockB::ReadTable(vector<vector<vector<vector<vector<double > > > > >* v, istream *table ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
 	for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
	  for(unsigned int i4=0;i4<v->at(i0)[i1][i2][i3].size();i4++){
	    *table >> v->at(i0)[i1][i2][i3][i4];
	    nn++;
	  }
	}
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int FastNLOBlockB::ReadTable(vector<vector<vector<vector<double > > > >* v, istream *table ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
 	for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
	  *table >> v->at(i0)[i1][i2][i3];
	  nn++;
	}
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int FastNLOBlockB::ReadTable(vector<vector<vector<double > > >* v, istream *table ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
	*table >> v->at(i0)[i1][i2];
	nn++;
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int FastNLOBlockB::ReadTable(vector<vector<double > >* v, istream *table ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      *table >> v->at(i0)[i1];
      nn++;
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int FastNLOBlockB::ReadTable(vector<double >* v, istream *table ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    *table >> v->at(i0);
    nn++;
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
void FastNLOBlockB::ResizeTable(vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int* dim5GetNxmaxFromDimI , int dim6 ){
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

void FastNLOBlockB::ResizeTable(vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5, int dim6 ){
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


void FastNLOBlockB::ResizeTable( vector<vector<vector<vector<vector<vector<double > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5 ){
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


void FastNLOBlockB::ResizeTable( vector<vector<vector<vector<vector<double > > > > >* v, int dim0 , int dim1, int dim2, int* dim3GetNxmaxFromDimI, int dim4 ){
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


void FastNLOBlockB::ResizeTable( vector<vector<vector<vector<vector<double > > > > >* v, int dim0 , int* dim1GetNxmaxFromDimI, int dim2, int dim3, int dim4 ){
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


void FastNLOBlockB::ResizeTable( vector<vector<vector<vector<vector<double > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4 ){
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

void FastNLOBlockB::ResizeTable( vector<vector<vector<vector<double > > > >* v, int dim0 , int dim1, int* dim2GetNxmaxFromDimI, int dim3 ){
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

void FastNLOBlockB::ResizeTable( vector<vector<vector<vector<double > > > >* v, int dim0 , int dim1, int dim2, int dim3 ){
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


void FastNLOBlockB::ResizeTable( vector<vector<vector<double > > >* v, int dim0 , int* dim1GetNxmaxFromDimI, int dim2 ){
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


void FastNLOBlockB::ResizeTable( vector<vector<vector<double > > >* v, int dim0 , int dim1, int dim2 ){
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


void FastNLOBlockB::ResizeTable( vector<vector<double > >*  v, int dim0 , int* dim1GetNxmaxFromDimI ){
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


void FastNLOBlockB::ResizeTable( vector<vector<double > >*  v, int dim0 , int dim1 ){
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


void FastNLOBlockB::ResizeTable( vector<double >* v, int dim0 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
  }
  else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //


int FastNLOBlockB::GetNxmax(int i){
   int nxmax = 0;
   switch (NPDFDim) {
   case 0: nxmax = Nxtot1[i];
      break;
   case 1: nxmax = ((int)pow((double)Nxtot1[i],2)+Nxtot1[i])/2;
      break;
   case 2: nxmax = Nxtot1[i]*Nxtot2[i];
      break;
   default: ;
   }
   return nxmax;
};


//________________________________________________________________________________________________________________ //


int FastNLOBlockB::GetTotalScalevars(){
   int totalscalevars=1;
   for(int scaledim=0;scaledim<NScaleDim;scaledim++){
      totalscalevars *= Nscalevar[scaledim];
   }
   return totalscalevars;
}

int FastNLOBlockB::GetTotalScalenodes(){
   int totalscalenodes=1;
   for(int scaledim=0;scaledim<NScaleDim;scaledim++){
      totalscalenodes *= Nscalenode[scaledim];
   }
   return totalscalenodes;
}


//________________________________________________________________________________________________________________ //


