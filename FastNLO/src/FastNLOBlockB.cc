// Author: Daniel Britzger
// DESY, 23/07/2011

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
  fname		= name;
  fNObsBins	= NObsBins;
}


FastNLOBlockB::FastNLOBlockB(const char* name , const int NObsBins , istream* table)
{
  fname		= name;
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
  *table >> NScaleDep;
  int NContrDescr;
  *table >> NContrDescr;
  // printf("# FastNLOBlockB::Read(): IDataFlag: %d, IAddMultFlag: %d, IContrFlag1: %d, IContrFlag2: %d, NScaleDep: %d\n",
  // 	 IDataFlag,IAddMultFlag,IContrFlag1,IContrFlag2,NScaleDep );
  CtrbDescript.resize(NContrDescr);
  char buffer[257];
  table->getline(buffer,256);
  for(int i=0;i<NContrDescr;i++){
    table->getline(buffer,256);
    CtrbDescript[i] = buffer;
    StripWhitespace(&CtrbDescript[i]);
  }

  int NCodeDescr;
  *table >> NCodeDescr;
  CodeDescript.resize(NCodeDescr);
  table->getline(buffer,256);
  for(int i=0;i<NCodeDescr;i++){
    table->getline(buffer,256);
    CodeDescript[i] = buffer;
    StripWhitespace(&CodeDescript[i]);
  }

  if(IDataFlag==1){
    *table >> Nuncorrel;
    UncDescr.resize(Nuncorrel);
    table->getline(buffer,256);
    for(int i=0;i<Nuncorrel;i++){
      table->getline(buffer,256);
      UncDescr[i] = buffer;
      StripWhitespace(&UncDescr[i]);
    }

    *table >> Ncorrel;
    CorDescr.resize(Ncorrel);
    table->getline(buffer,256);
    for(int i=0;i<Ncorrel;i++){
      table->getline(buffer,256);
      CorDescr[i] = buffer;
      StripWhitespace(&CorDescr[i]);
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
      StripWhitespace(&UncDescr[i]);
    }
    *table >> Ncorrel;
    CorDescr.resize(Ncorrel);
    table->getline(buffer,256);
    for(int i=0;i<Ncorrel;i++){
      table->getline(buffer,256);
      CorDescr[i] = buffer;
      StripWhitespace(&CorDescr[i]);
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

    int NscaleDescript;
    ScaleDescript.resize(NScaleDim);
    for(int i=0;i<NScaleDim;i++){
      *table >> NscaleDescript;
      ScaleDescript[i].resize(NscaleDescript);
      table->getline(buffer,256);
      for(int j=0;j<NscaleDescript;j++){
	table->getline(buffer,256);
	ScaleDescript[i][j] = buffer;
	StripWhitespace(&ScaleDescript[i][j]);
      }
    }


    if ( NScaleDep < 3 ) {
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

    }
      
    if ( NScaleDep >= 3 ) {
      int nn3 = 0;

      nn3 += ReadFlexibleVector  ( &ScaleNodeScale1 , table );
      nn3 += ReadFlexibleVector  ( &ScaleNodeScale2 , table );

      nn3 += ReadFlexibleVector  ( &SigmaTildeMuIndep , table , true );
      if ( NScaleDep==3 || NScaleDep==5) {
	 nn3 += ReadFlexibleVector  ( &SigmaTildeMuFDep , table , true );
	 nn3 += ReadFlexibleVector  ( &SigmaTildeMuRDep , table , true );
      }
      nn3 += ReadFlexibleVector  ( &SigmaRefMixed , table , true );
      nn3 += ReadFlexibleVector  ( &SigmaRef_s1 , table , true );
      nn3 += ReadFlexibleVector  ( &SigmaRef_s2 , table , true );
         
      ResizeFlexibleVector( &PdfLcMuVar  , &SigmaTildeMuIndep );
      AlphasTwoPi.resize(ScaleNodeScale1.size());
      for ( unsigned int i=0; i<AlphasTwoPi.size() ; i++ ){
	AlphasTwoPi[i].resize(ScaleNodeScale1[i].size());
	for ( unsigned int j=0; j<AlphasTwoPi[i].size() ; j++ ){
	  AlphasTwoPi[i][j].resize(ScaleNodeScale2[i].size());
	}
      }

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



void FastNLOBlockB::Print(const int ic, const int iprint){
  if (iprint==0) {
    printf(" # Contribution %1i:\n",ic);
    for(unsigned int i=0;i<CtrbDescript.size();i++){
      printf(" #   %s\n",CtrbDescript[i].data());
    }
    printf(" #   No. of events: %16llu\n",Nevt);
    printf(" #   provided by:\n");
    for(unsigned int i=0;i<CodeDescript.size();i++){
      printf(" #   %s\n",CodeDescript[i].data());
    }
    if ( NScaleDep<3 ) {
      printf(" #   Scale dimensions: %1i\n",NScaleDim);
      for(int i=0;i<NScaleDim;i++){
	for(unsigned int j=0;j<ScaleDescript[i].size();j++){
	  printf(" #     Scale description for dimension %1i:          %s\n",i+1,ScaleDescript[i][j].data());
	}
	printf(" #     Number of scale variations for dimension %1i: %1i\n",i+1,Nscalevar[i]);
	printf(" #     Available scale settings for dimension %1i:\n",i+1);
	for(int k=0;k<Nscalevar[i];k++){
	  printf(" #       Scale factor number %1i:                   % #10.4f\n",k+1,ScaleFac[i][k]);
	}
	printf(" #     Number of scale nodes for dimension %1i:      %1i\n",i+1,Nscalenode[i]);
      }
    }
  } else {
    if (ic==1) {
      printf("\n *****************************************\n");
      printf(" * fastNLO Table: Block B\n");
      printf(" *****************************************\n");
    }
    //  printf("  B0  fNNObsBins                        %10i\n",fNObsBins);
    printf("  B0  ISep                              %10i\n",tablemagicno);
    printf("  B0    IXsectUnits(%1i)                  %10i\n",ic,IXsectUnits);
    printf("  B0    IDataFlag(%1i)                    %10i\n",ic,IDataFlag);
    printf("  B0    IAddMultFlag(%1i)                 %10i\n",ic,IAddMultFlag);
    printf("  B0    IContrFlag1(%1i)                  %10i\n",ic,IContrFlag1);
    printf("  B0    IContrFlag2(%1i)                  %10i\n",ic,IContrFlag2);
    printf("  B0    NScaleDep(%1i)                    %10i\n",ic,NScaleDep);
    printf("  B0    NContrDescr(%1i)                  %10zi\n",ic,CtrbDescript.size());
    for(unsigned int i=0;i<CtrbDescript.size();i++){
      printf("  B0      CtrbDescript(%1i,%1i)             %s\n",ic,i+1,CtrbDescript[i].data());
    }
    printf("  B0    NCodeDescr(%1i)                   %10zi\n",ic,CodeDescript.size());
    for(unsigned int i=0;i<CodeDescript.size();i++){
      printf("  B0      CodeDescript(%1i,%1i)             %s\n",ic,i+1,CodeDescript[i].data());
    }

    if(IDataFlag==1){
      printf("  B0  Nuncorrel                         %10i\n",Nuncorrel);
      printf("  B0  Ncorrel                           %10i\n",Ncorrel);
      printf("  B0  NErrMatrix                        %10i\n",NErrMatrix);
      printf("  B0  some more output could be printed here (IDataFlag==1).\n");
    }
    if(IAddMultFlag==1){
      printf("  B0  some more output could be printed here (IAddMultFlag==1).\n");
    }

    if(!(IDataFlag==1) && !(IAddMultFlag==1)){ // that's the usual case
      printf("  B0    IRef(%1i)                         %10i\n",ic,IRef);
      printf("  B0    IScaleDep(%1i)                    %10i\n",ic,IScaleDep);
      printf("  B0    Nevt(%1i)                         %16llu\n",ic,Nevt);
      printf("  B0    Npow(%1i)                         %10i\n",ic,Npow);
      printf("  B0    NPDF(%1i)                         %10i\n",ic,NPDF);
      if(NPDF>0){
	for(int i=0;i<NPDF;i++){
	  printf("  B0      NPDFPDG(%1i,%1i)                  %10i\n",ic,i+1,NPDFPDG[i]);
	}
      }
      printf("  B0    NPDFDim(%1i)                      %10i\n",ic,NPDFDim);
      printf("  B0    NFragFunc(%1i)                    %10i\n",ic,NFragFunc);
      if(NFragFunc>0){
	for(int i=0;i<NFragFunc;i++){
	  printf("  B0      NFFPDG(%1i,%1i)                   %10i\n",ic,i+1,NFFPDG[i]);
	}
      }
      printf("  B0    NFFDim(%1i)                       %10i\n",ic,NFFDim);
      printf("  B0    NSubproc(%1i)                     %10i\n",ic,NSubproc);
      printf("  B0    IPDFdef1(%1i)                     %10i\n",ic,IPDFdef1);
      printf("  B0    IPDFdef2(%1i)                     %10i\n",ic,IPDFdef2);
      printf("  B0    IPDFdef3(%1i)                     %10i\n",ic,IPDFdef3);
      printf("  B0    NScales(%1i)                      %10i\n",ic,NScales);
      printf("  B0    NScaleDim(%1i)                    %10i\n",ic,NScaleDim);
      for(int i=0;i<NScales;i++){
	printf("  B0      NScales(%1i,%1i)                  %10i\n",ic,i+1,Iscale[i]);
      }
      if ( NScaleDep<3 ) {
	for(int i=0;i<NScaleDim;i++){
	  printf("  B0      NScaleDescript(%1i,%1i)           %10i\n",ic,i+1,NScaleDim);
	  for(unsigned int j=0;j<ScaleDescript[i].size();j++){
	    printf("  B0        ScaleDescript(%1i,%1i,%1i)        %s\n",ic,i+1,j+1,ScaleDescript[i][j].data());
	  }
	  printf("  B0      NScaleVar(%1i,%1i)                %10i\n",ic,i+1,Nscalevar[i]);
	  printf("  B0      NScaleNode(%1i,%1i)               %10i\n",ic,i+1,Nscalenode[i]);
	}
      }
      if (ic>1) {
	printf(" #########################################\n");
      }
    }
  }
}



//______________________________________________________________________________
/*
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
*/
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
/*
//________________________________________________________________________________________________________________ //
int FastNLOBlockB::ReadFlexibleVector(vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, istream *table , bool nProcLast ){
  int nn = 0;
  int size = 0;
  *table >> size; nn++;
  v->resize(size);
  for(unsigned int i0=0;i0<v->size();i0++){
    nn += ReadFlexibleVector(&(v->at(i0)),table,nProcLast);
  }
  return nn;
}
//________________________________________________________________________________________________________________ //
int FastNLOBlockB::ReadFlexibleVector(vector<vector<vector<vector<vector<vector<double > > > > > >* v, istream *table , bool nProcLast ){
  int nn = 0;
  int size = 0;
  *table >> size; nn++;
  v->resize(size);
  for(unsigned int i0=0;i0<v->size();i0++){
    nn += ReadFlexibleVector(&(v->at(i0)),table,nProcLast);
  }
  return nn;
}
//________________________________________________________________________________________________________________ //
int FastNLOBlockB::ReadFlexibleVector(vector<vector<vector<vector<vector<double > > > > >* v, istream *table , bool nProcLast ){
  int nn = 0;
  int size = 0;
  *table >> size; nn++;
  v->resize(size);
  for(unsigned int i0=0;i0<v->size();i0++){
    nn += ReadFlexibleVector(&(v->at(i0)),table,nProcLast);
  }
  return nn;
}
//________________________________________________________________________________________________________________ //
int FastNLOBlockB::ReadFlexibleVector(vector<vector<vector<vector<double > > > >* v, istream *table , bool nProcLast ){
  int nn = 0;
  int size = 0;
  *table >> size; nn++;
  v->resize(size);
  for(unsigned int i0=0;i0<v->size();i0++){
    nn += ReadFlexibleVector(&(v->at(i0)),table,nProcLast);
  }
  return nn;
}
//________________________________________________________________________________________________________________ //
int FastNLOBlockB::ReadFlexibleVector(vector<vector<vector<double > > >* v, istream *table , bool nProcLast ){
  int nn = 0;
  int size = 0;
  *table >> size; nn++;
  v->resize(size);
  for(unsigned int i0=0;i0<v->size();i0++){
    nn += ReadFlexibleVector(&(v->at(i0)),table,nProcLast);
  }
  return nn;
}
//________________________________________________________________________________________________________________ //
int FastNLOBlockB::ReadFlexibleVector(vector<vector<double > >* v, istream *table , bool nProcLast ){
  int nn = 0;
  int size = 0;
  *table >> size; nn++;
  v->resize(size);
  for(unsigned int i0=0;i0<v->size();i0++){
    nn += ReadFlexibleVector(&(v->at(i0)),table,nProcLast);
  }
  return nn;
}
*/
//________________________________________________________________________________________________________________ //
int FastNLOBlockB::ReadFlexibleVector(vector<double >* v, istream *table , bool nProcLast ){
  int nn = 0;
  if ( !nProcLast ){
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

// //________________________________________________________________________________________________________________ //
// void FastNLOBlockB::ResizeFlexibleVector(vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, vector<vector<vector<vector<vector<vector<vector<double > > > > > > >*nom ){
//   // resize vector v to size of vector nom
//   v->resize(nom->size());
//   for ( unsigned int i = 0 ; i<v->size() ; i++ ){
//     ResizeFlexibleVector(&((*v)[i]),&((*nom)[i]));
//   }
// }
// //________________________________________________________________________________________________________________ //
// void FastNLOBlockB::ResizeFlexibleVector(vector<vector<vector<vector<vector<vector<double > > > > > >* v, vector<vector<vector<vector<vector<vector<double > > > > > >*nom ){
//   // resize vector v to size of vector nom
//   v->resize(nom->size());
//   for ( unsigned int i = 0 ; i<v->size() ; i++ ){
//     ResizeFlexibleVector(&((*v)[i]),&((*nom)[i]));
//   }
// }
// //________________________________________________________________________________________________________________ //
// void FastNLOBlockB::ResizeFlexibleVector(vector<vector<vector<vector<vector<double > > > > >* v, vector<vector<vector<vector<vector<double > > > > >*nom ){
//   // resize vector v to size of vector nom
//   v->resize(nom->size());
//   for ( unsigned int i = 0 ; i<v->size() ; i++ ){
//     ResizeFlexibleVector(&((*v)[i]),&((*nom)[i]));
//   }
// }
// //________________________________________________________________________________________________________________ //
// void FastNLOBlockB::ResizeFlexibleVector(vector<vector<vector<vector<double > > > >* v, vector<vector<vector<vector<double > > > >*nom ){
//   // resize vector v to size of vector nom
//   v->resize(nom->size());
//   for ( unsigned int i = 0 ; i<v->size() ; i++ ){
//     ResizeFlexibleVector(&((*v)[i]),&((*nom)[i]));
//   }
// }
// //________________________________________________________________________________________________________________ //
// void FastNLOBlockB::ResizeFlexibleVector(vector<vector<vector<double > > >* v, vector<vector<vector<double > > >*nom ){
//   // resize vector v to size of vector nom
//   v->resize(nom->size());
//   for ( unsigned int i = 0 ; i<v->size() ; i++ ){
//     ResizeFlexibleVector(&((*v)[i]),&((*nom)[i]));
//   }
// }
// //________________________________________________________________________________________________________________ //
// void FastNLOBlockB::ResizeFlexibleVector(vector<vector<double > >* v, vector<vector<double > >*nom ){
//   // resize vector v to size of vector nom
//   v->resize(nom->size());
//   for ( unsigned int i = 0 ; i<v->size() ; i++ ){
//     ResizeFlexibleVector(&((*v)[i]),&((*nom)[i]));
//   }
// }
// //________________________________________________________________________________________________________________ //
// void FastNLOBlockB::ResizeFlexibleVector(vector<double >* v, vector<double >*nom ){
//   // resize vector v to size of vector nom
//   v->resize(nom->size());
// }


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

void FastNLOBlockB::StripWhitespace(string* s){
  string fastlast = &(*s)[s->size()-1];
  while ( !fastlast.compare(" ")){
    string::iterator it = s->end();
    s->erase(it-1);
    fastlast = &(*s)[s->size()-1];
  }
}

//________________________________________________________________________________________________________________ //
