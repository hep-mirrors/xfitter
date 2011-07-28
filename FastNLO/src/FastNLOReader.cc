// Author: Daniel Britzger
// DESY, 23/07/2011

//  Version 0.1, 

////////////////////////////////////////////////////////////////////////
//
// FastNLOReader
//
//  ... more documentation ...
//
//

#include "FastNLOReader.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include "LHAPDF.h"

using namespace std;

//ClassImp(FastNLOReader)

//______________________________________________________________________________



FastNLOReader::FastNLOReader(void)
{
  //
  // do not call the standard constructor
  // 
  cout << "FastNLOReader::FastNLOReader. Please set filename and run ReadTable()! "<<endl;
  fScalevar	= 0;
  SetPDFInterface(FastNLOReader::kLHAPDF);
  SetAlphasEvolution(FastNLOReader::kGRV);
}




FastNLOReader::FastNLOReader(string filename)
{

  printf(" ***************************************************************** \n");
  printf(" *  \n");
  printf(" *  FastNLO Reader - version 0.1\n");
  printf(" *  \n");
  printf(" *  This code is a reader for FastNLO tables, that were\n");
  printf(" *  calculated using nlojet++ v4.1.3 and FastNLO v2.0 (or higher).\n"); 
  printf(" *  \n"); 
  printf(" *  FastNLO was developed by\n"); 
  printf(" *    Thomas Kluge, Klaus Rabbertz, Markus Wobisch\n");
  printf(" *    (publication in preparation)\n");
  printf(" *  \n");
  printf(" *  Modifications to FastNLO and developement and maintenance\n");
  printf(" *  of this code by\n");
  printf(" *    Daniel Britzger\n");
  printf(" *    daniel.britzger@desy.de \n");
  printf(" *  \n");
  printf(" *  for NLOJET++ please cite \n");
  printf(" *    Z. Nagy, Phys. Rev. Lett. 88, 122003 (2002),\n");
  printf(" *    Z. Nagy, Phys. Rev. Lett. 88, 122003 (2002), \n");
  printf(" *  \n");
  printf(" ***************************************************************** \n");

  SetFilename(filename);
  fScalevar	= 0;
  ReadTable();
  SetPDFInterface(FastNLOReader::kLHAPDF);
  SetAlphasEvolution(FastNLOReader::kGRV);
  
  //Print();
}




//______________________________________________________________________________



double FastNLOReader::SetScaleVariation(int scalevar){
  
  // ------------------------------------------------
  //   Set the scalevariation factor for detemining the
  //   'theory'-error. Usually, you have tables stored with
  //   factors of 0.5, 1 and 2 times the nominal scale.
  //     corresponding to:
  //     scalevar -> scalefactor
  //        '0'   ->   1
  //        '1'   ->   0.5
  //        '2'   ->   2
  //   This method returns the scalefactor correspoding to
  //   the choosen 'scalevar'.
  // ------------------------------------------------

  if ( BlockB_NLO->NScaleDep == 3 ){
    printf("FastNLOReader::SetScaleVariation(). Info: This is a MuVar table, therefore, you can choose all possible scale variations. Your Scalevar has to be '0'.\n");
  }
  

  if (  scalevar >= BlockB_NLO->Nscalevar[0]  ){
    printf("Warning in FastNLOReader::SetScaleVariation. This table has only %d scalevariations stored. You wanted to acces number %d. Using '0' instead.\n", BlockB_NLO->Nscalevar[0] ,scalevar );
    fScalevar	= 0;
    return BlockB_NLO->ScaleFac[0][0] ;
  }
  
  fScalevar	= scalevar;
  printf(" * FastNLOReader::SetScaleVariation. Scalefactor of %4.2f for the nominal scale is chosen.\n",BlockB_NLO->ScaleFac[0][fScalevar]);
  return BlockB_NLO->ScaleFac[0][fScalevar];
  
}




//______________________________________________________________________________



void FastNLOReader::ReadTable(void)
{
  // ---- check if file exists ----- //
  FILE* fp = fopen(ffilename.c_str(), "r");
  if (fp) {
    fclose(fp);
  } else {
    printf("Error. File does not exists. Was looking for: %s. Exiting.\n",ffilename.c_str());
    exit(1);
  } 
  
  // open stream
  ifstream* instream = new ifstream(ffilename.c_str(),ios::in);

  ReadBlockA1(instream);
  ReadBlockA2(instream);
  //fnloBlockA1 *blocka1 = GetBlockA1();  
  //int nblocks = blocka1->GetNcontrib()+blocka1->GetNdata();
  int nblocks	= Ncontrib + Ndata;
  printf(" * FastNLOReader::ReadTable(). Reading %d B-Blocks.\n",nblocks);
  for(int i=0;i<nblocks;i++){
    FastNLOBlockB* blockb	= new FastNLOBlockB( "ReadingBlockB", NObsBin , instream );
    if ( blockb->IContrFlag1==1 && blockb->IContrFlag2==1 && blockb->IRef==0 ){
      if ( blockb->NScaleDep == 0 )   blockb->SetName("BlockB. LO. v2.0.");
      if ( blockb->NScaleDep == 2 )   blockb->SetName("BlockB. LO. v2.0-TwoScales.");
      if ( blockb->NScaleDep == 3 )   blockb->SetName("BlockB. LO. v2.0 MuVar.");
      BlockB_LO		= blockb;
    }
    else if ( blockb->IContrFlag1==1 && blockb->IContrFlag2==2 && blockb->IRef==0 ){
      if ( blockb->NScaleDep == 0 )   blockb->SetName("BlockB. NLO. v2.0.");
      if ( blockb->NScaleDep == 2 )   blockb->SetName("BlockB. NLO. v2.0-TwoScales.");
      if ( blockb->NScaleDep == 3 )   blockb->SetName("BlockB. NLO. v2.0 MuVar.");
      BlockB_NLO	= blockb;
    }
    else if ( blockb->IContrFlag1==1 && blockb->IContrFlag2==1 && blockb->IRef==1 ){
      if ( blockb->NScaleDep == 0 )   blockb->SetName("BlockB. LO Reference. v2.0.");
      if ( blockb->NScaleDep == 2 )   blockb->SetName("BlockB. LO Reference. v2.0-TwoScales.");
      if ( blockb->NScaleDep == 3 )   blockb->SetName("BlockB. LO Reference. v2.0 MuVar.");
      BlockB_LO_Ref		= blockb;
    }
    else if ( blockb->IContrFlag1==1 && blockb->IContrFlag2==2 && blockb->IRef==1 ){
      if ( blockb->NScaleDep == 0 )   blockb->SetName("BlockB. NLO Reference. v2.0.");
      if ( blockb->NScaleDep == 2 )   blockb->SetName("BlockB. NLO Reference. v2.0-TwoScales.");
      if ( blockb->NScaleDep == 3 )   blockb->SetName("BlockB. NLO Reference. v2.0 MuVar.");
      BlockB_NLO_Ref	= blockb;
    }
    else {
      printf("FastNLOReader::ReadTable(). Error in initializing the 'B'-Blocks. Exiting.\n");
      exit(1);      
    }
  }

  if ( (BlockB_LO_Ref == NULL || BlockB_NLO_Ref == NULL) && BlockB_LO->NScaleDep==0 ){
    printf("No Reference Tables were found.\n");
  }
  
  if ( BlockB_LO == NULL ){
    printf("ERROR. Could not find any LO Calculation (BlockB_LO).\n");exit(1);
  }
  
  //NPDFDim	= BlockB_LO->NPDFDim;

}

//______________________________________________________________________________


void FastNLOReader::ReadBlockA1(istream *table){
   table->peek();
   if (table->eof()){
      printf("BlockA1::Read: Cannot read from file.\n");
      return;
   }

   int key = 0;
   *table >> key;
   if(key != tablemagicno){
      printf("fnloBlockA1::Read: At beginning of block found %d instead of %d.\n",key,tablemagicno);
      return;
   };
   *table >> Itabversion;
   *table >> ScenName;
   *table >> Ncontrib;
   *table >> Nmult;
   *table >> Ndata;
   *table >> NuserString;
   *table >> NuserInt;
   *table >> NuserFloat;
   *table >> Imachine;
   key=0;
   *table >> key;
   if(key != tablemagicno){
      printf("fnloBlockA1::Read: At end of block found %d instead of %d.\n",key,tablemagicno);
      return;
   };
   // Put magic number back
   for(int i=0;i<(int)(log10((double)key)+1);i++){
      table->unget();
   }
}


//______________________________________________________________________________


void FastNLOReader::ReadBlockA2(istream *table){
   table->peek();
   if (table->eof()){
      printf("fnloBlockA2::Read: Cannot read from file.\n");
      return;
   }

   int key = 0;
   *table >> key;
   if(key != tablemagicno){
      printf("fnloBlockA2::Read: At beginning of block found %d instead of %d.\n",key,tablemagicno);
      return;
   };

   *table >> Ipublunits;
   *table >> NScDescript;
   ScDescript.resize(NScDescript);
   char buffer[257];
   table->getline(buffer,256);
   for(int i=0;i<NScDescript;i++){
      table->getline(buffer,256);
      ScDescript[i] = buffer;
      //      StripWhitespace(ScDescript[i]);
   }

   *table >> Ecms;
   *table >> ILOord;
   *table >> NObsBin;
   *table >> NDim;
   DimLabel.resize(NDim);
   table->getline(buffer,256);
   for(int i=0;i<NDim;i++){
      table->getline(buffer,256);
      DimLabel[i] = buffer;
      //      StripWhitespace(DimLabel[i]);
   }

   IDiffBin.resize(NDim);
   for(int i=0;i<NDim;i++){
      *table >>  IDiffBin[i];
   }
   LoBin.resize(NObsBin);
   UpBin.resize(NObsBin);
   //KR: Set rapidity index also when reading a table
   RapIndex.push_back(0);
   //   int irap = 0;
   for(int i=0;i<NObsBin;i++){
      LoBin[i].resize(NDim);
      UpBin[i].resize(NDim);
      for(int j=0;j<NDim;j++){
         *table >>  LoBin[i][j];
         if(IDiffBin[j]==2) *table >>  UpBin[i][j];
      }
      //      cout << "iobs1: " << i << ", LoBin i: " << LoBin[i][1] << endl;
      if ( i > 0 ) {
        if ( LoBin[i][1] != LoBin[i-1][1] ) {
          //      cout << "iobs2: " << i << ", LoBin i-1: " << LoBin[i-1][1] << ", LoBin i: " << LoBin[i][1] << endl;
          RapIndex.push_back(i);
          //      irap++;
          //      cout << "irap: " << irap << ", RapIndex: " << RapIndex[irap] << endl;
        }
      }
   }

   BinSize.resize(NObsBin);
   for(int i=0;i<NObsBin;i++){
     *table >> BinSize[i];
   }

   *table >> INormFlag;
   if(INormFlag>1){
      *table >> DenomTable;
   }
   if(INormFlag>0){
      IDivLoPointer.resize(NObsBin);
      IDivUpPointer.resize(NObsBin);
      for(int i=0;i<NObsBin;i++){
         *table >> IDivLoPointer[i];
         *table >> IDivUpPointer[i];
      }
   }

   key=0;
   *table >> key;
   if(key != tablemagicno){
      printf("fnloBlockA2::Read: At end of block found %d instead of %d.\n",key,tablemagicno);
      return;
   };
   // Put magic number back
   for(int i=0;i<(int)(log10((double)key)+1);i++){
      table->unget();
   }
}



//______________________________________________________________________________


void FastNLOReader::Print(){
  PrintBlockA1();
  PrintBlockA2();
  if ( BlockB_LO      ) BlockB_LO	->Print();
  if ( BlockB_NLO     ) BlockB_NLO	->Print();
  if ( BlockB_LO_Ref  ) BlockB_LO_Ref	->Print();
  if ( BlockB_NLO_Ref ) BlockB_NLO_Ref	->Print();
}



//______________________________________________________________________________


void FastNLOReader::PrintBlockA1(){
  printf("\n **************** FastNLO Table: BlockA1 ****************\n\n");
  printf(" A1  tablemagicno                  %d\n",tablemagicno);
  printf(" A1  Itabversion                   %d\n",Itabversion);
  printf(" A1  ScenName                       %s\n",ScenName.data());
  printf(" A1  Ncontrib                      %d\n",Ncontrib);
  printf(" A1  Nmult                         %d\n",Nmult);
  printf(" A1  Ndata                         %d\n",Ndata);
  printf(" A1  NuserString                   %d\n",NuserString);
  printf(" A1  NuserInt                      %d\n",NuserInt);
  printf(" A1  NuserFloat                    %d\n",NuserFloat);
  printf(" A1  Imachine                      %d\n",Imachine);
  printf("\n ********************************************************\n\n");
}



//______________________________________________________________________________


void FastNLOReader::PrintBlockA2(){
  printf("\n **************** FastNLO Table: BlockA2 ****************\n\n");
  printf(" A2  Ipublunits                    %d\n",Ipublunits);
  printf(" A2  NScDescript                   %d\n",NScDescript);
  for(int i=0;i<NScDescript;i++){
    printf(" A2  ScDescript[%d]                 %s\n",i,ScDescript[i].data());
  }
  printf(" A2  Ecms                          %7.4f\n",Ecms);
  printf(" A2  ILOord                        %d\n",ILOord);
  printf(" A2  NDim                          %d\n",NDim);
  for(int i=0;i<NDim;i++){
    printf(" A2   - DimLabel[%d]                %s\n",i,DimLabel[i].data());
  }
  for(int i=0;i<NDim;i++){
    printf(" A2   - IDiffBin[%d]               %d\n",i,IDiffBin[i]);
  }
  printf(" A2  NObsBin                       %d\n",NObsBin);
  for(int i=0;i<NObsBin;i++){
    for(int j=0;j<NDim;j++){
      printf(" A2   -  - LoBin[%d][%d]             %7.4f\n", i,j,LoBin[i][j]);
      if(IDiffBin[j]==2)
        printf(" A2   -  - UpBin[%d][%d]             %7.4f\n", i,j,UpBin[i][j]);
    }
   }
   for(int i=0;i<NObsBin;i++){
     printf(" A2   - BinSize[%d]                %7.4f\n", i,BinSize[i]);
   }
   printf(" A2  INormFlag                     %d\n",INormFlag);

   if(INormFlag>1){
     printf(" A2  DenomTable                    %s\n",DenomTable.data());
   }
   if(INormFlag>0){
      for(int i=0;i<NObsBin;i++){
        printf(" A2   - IDivLoPointer[%d]               %d\n",i,IDivLoPointer[i]);
        printf(" A2   - IDivUpPointer[%d]               %d\n",i,IDivUpPointer[i]);
      }
   }
   printf("\n ********************************************************\n\n");
}


//______________________________________________________________________________


vector < double > FastNLOReader::GetXSection( ){
  // check, if x-section is already calculated...
  
  if ( XSection.empty() ){
    CalcCrossSection();
  }
  
  if ( BlockB_NLO->NScaleDep == 3 )
    return XSectionMuVar;
  else return XSection;
}


//______________________________________________________________________________


vector < double > FastNLOReader::GetReferenceXSection( ){
  // check, if x-section is already calculated...
  
  if ( XSectionRef.empty() ){
    CalcReferenceCrossSection();
  }
  
  if ( BlockB_NLO->NScaleDep == 3 )
    return XSectionRefQ2;
  else return XSectionRef;
}


//______________________________________________________________________________


void FastNLOReader::CalcReferenceCrossSection( ){

  XSectionRef.clear();
  XSection2ScalesRef.clear();
  XSectionRef.resize(NObsBin);
  XSection2ScalesRef.resize(NObsBin);
  //XSectionMuVarRef.clear();
  //XSectionMuVarRef.resize(NObsBin);

  XSectionRefMixed.clear();		  
  XSectionRefQ2.clear();
  XSectionRefMufQ2MuRMixed.clear();
  XSectionRefMixed.resize(NObsBin);		  
  XSectionRefQ2.resize(NObsBin);
  XSectionRefMufQ2MuRMixed.resize(NObsBin);

  if ( BlockB_LO_Ref && BlockB_NLO_Ref ){

    for(int i=0;i<NObsBin;i++){
      for(int l=0;l<BlockB_LO_Ref->NSubproc;l++){ 
	XSectionRef[i] +=  BlockB_LO_Ref->SigmaTilde[i][fScalevar][0][0][l] * BinSize[i];
      }
      for(int l=0;l<BlockB_NLO_Ref->NSubproc;l++){ 
	XSectionRef[i] +=  BlockB_NLO_Ref->SigmaTilde[i][fScalevar][0][0][l] * BinSize[i];
      }
    }
    
  }
  
  if ( BlockB_LO->NScaleDep == 2 ){
    // todo
  }
  
  if ( BlockB_LO->NScaleDep == 3 ){
    for(int i=0;i<NObsBin;i++){
      for(int n=0;n<BlockB_NLO->NSubproc;n++) {
	XSectionRefMixed[i]		+= BlockB_LO ->SigmaRefMixed[i][n] * BinSize[i];
	XSectionRefQ2[i]		+= BlockB_LO ->SigmaRefQ2[i][n] * BinSize[i];
	XSectionRefMufQ2MuRMixed[i]	+= BlockB_LO ->SigmaRefMufQ2MuRMixed[i][n] * BinSize[i];
      }
      for(int n=0;n<BlockB_NLO->NSubproc;n++) {
	XSectionRefMixed[i]		+= BlockB_NLO->SigmaRefMixed[i][n] * BinSize[i];
	XSectionRefQ2[i]		+= BlockB_NLO->SigmaRefQ2[i][n] * BinSize[i];
	XSectionRefMufQ2MuRMixed[i]	+= BlockB_NLO->SigmaRefMufQ2MuRMixed[i][n] * BinSize[i];
      }
    }
  }

}


//______________________________________________________________________________


void FastNLOReader::CalcCrossSection( ){
  
  // todo! Check if pdflcs and alphas cache is filled...

  XSection_LO.clear();
  XSection.clear();
  XSection2Scales_LO.clear();
  XSection2Scales.clear();
  XSectionMuVar_LO.clear();
  XSectionMuVar.clear();

  XSection_LO.resize(NObsBin);
  XSection.resize(NObsBin);
  XSection2Scales_LO.resize(NObsBin);
  XSection2Scales.resize(NObsBin);
  XSectionMuVar_LO.resize(NObsBin);
  XSectionMuVar.resize(NObsBin);
  
  for(int i=0;i<NObsBin;i++){
    int nxmax = BlockB_LO->GetNxmax(i);
    
    for(int j=0;j<BlockB_LO->GetTotalScalenodes();j++){
      int scalenode1 = j;
      int scalenode2 = j;
      if (BlockB_LO->NScaleDim>1){          
	scalenode1 = j / BlockB_LO->Nscalenode[1];
	scalenode2 = j % BlockB_LO->Nscalenode[1];
      }
      
      for(int k=0;k<nxmax;k++){ 
	// LO Block
	for(int l=0;l<BlockB_LO->NSubproc;l++){ 
	  XSection[i]		+=  BlockB_LO->SigmaTilde[i][fScalevar][j][k][l] *  BlockB_LO->AlphasTwoPi_v20[i][scalenode2]  *  BlockB_LO->PdfLc[i][scalenode2][k][l] * BinSize[i];
	  XSection_LO[i]	+=  BlockB_LO->SigmaTilde[i][fScalevar][j][k][l] *  BlockB_LO->AlphasTwoPi_v20[i][scalenode2]  *  BlockB_LO->PdfLc[i][scalenode2][k][l] * BinSize[i];
	  //printf("%15.13f     %9.5f         %16.14f\n",BlockB_LO->SigmaTilde[i][fScalevar][j][k][l] ,  BlockB_LO->AlphasTwoPi_v20[i][scalenode2]  ,  BlockB_LO->PdfLc[i][scalenode2][k][l]);
	}
	// NLO Block
	for(int l=0;l<BlockB_NLO->NSubproc;l++){ 
	  XSection[i]		+=  BlockB_NLO->SigmaTilde[i][fScalevar][j][k][l] *  BlockB_NLO->AlphasTwoPi_v20[i][scalenode2]  *  BlockB_NLO->PdfLc[i][scalenode2][k][l] * BinSize[i];
	}
      }
    }


    if ( BlockB_LO->NScaleDep == 2 ){
      for(int jNodeR=0;jNodeR<BlockB_LO->NscalenodeScale1;jNodeR++){
	for(int iNodeF=0;iNodeF<BlockB_LO->NscalenodeScale2;iNodeF++){
	  for(int x=0;x<nxmax;x++){ 
	    // LO Block
	    for(int n=0;n<BlockB_LO->NSubproc;n++){ 
	      XSection2Scales[i]	+=  BlockB_LO->SigmaTilde2Scales[i][fScalevar][fScalevar][jNodeR][iNodeF][x][n] * BlockB_LO->AlphasTwoPi2Scales[i][jNodeR]   *  BlockB_LO->PdfLc2Scales[i][iNodeF][x][n] * BinSize[i];
	      XSection2Scales_LO[i]	+=  BlockB_LO->SigmaTilde2Scales[i][fScalevar][fScalevar][jNodeR][iNodeF][x][n] * BlockB_LO->AlphasTwoPi2Scales[i][jNodeR]   *  BlockB_LO->PdfLc2Scales[i][iNodeF][x][n] * BinSize[i];
	    }
	    // NLO Block
	    for(int n=0;n<BlockB_NLO->NSubproc;n++){ 
	      XSection2Scales[i]	+=  BlockB_NLO->SigmaTilde2Scales[i][fScalevar][fScalevar][jNodeR][iNodeF][x][n] * BlockB_NLO->AlphasTwoPi2Scales[i][jNodeR]   *  BlockB_NLO->PdfLc2Scales[i][iNodeF][x][n] * BinSize[i];
	    }
	  }
	}
      }
   }

    if ( BlockB_LO->NScaleDep == 3 ){
      
      for(int jQ=0;jQ<BlockB_LO->NscalenodeScaleQ;jQ++){
	for(int jPt=0;jPt<BlockB_LO->NscalenodeScalePt;jPt++){
		    double Q2   = BlockB_LO->ScaleNodeQ[i][jQ]*BlockB_LO->ScaleNodeQ[i][jQ];
		    double Pt   = BlockB_LO->ScaleNodePt[i][jPt];
		    
		    //
		    //   TODO
		    //
		    // todo! here the formula for mu_f and mu_r calculation has to be implemented
		    // 		    double mur2 = (ScaleNodeQ[i][jQ]**2 + ScaleNodePt[i][jPt]*ScaleNodePt[i][jPt] )/ 2.;
		    // 		    double muf2 = (ScaleNodeQ[i][jQ]**2 + ScaleNodePt[i][jPt]*ScaleNodePt[i][jPt] )/ 2.;
		    
    		    double mur2 = Q2;
    		    double muf2 = Q2;

 		    for(int x=0;x<nxmax;x++){ 
		      
		      //  -----  DIS ---- //
		      if ( BlockB_LO->NPDFDim == 0 ) {
			// LO Block
			for(int n=0;n<BlockB_LO->NSubproc;n++){ 
			  double as	= BlockB_LO->AlphasTwoPi[i][jQ][jPt];
			  double pdflc	= BlockB_LO->PdfLcMuVar[i][x][jQ][jPt][n];
			  XSectionMuVar[i]	+=  BlockB_LO->SigmaTildeMuIndep[i][x][jQ][jPt][n] *                     as * pdflc * BinSize[i];
			  XSectionMuVar[i]	+=  BlockB_LO->SigmaTildeMuFDep [i][x][jQ][jPt][n] * std::log(muf2/Q2) * as * pdflc * BinSize[i];
			  XSectionMuVar[i]	+=  BlockB_LO->SigmaTildeMuRDep [i][x][jQ][jPt][n] * std::log(mur2/Q2) * as * pdflc * BinSize[i];

			  XSectionMuVar_LO[i]	+=  BlockB_LO->SigmaTildeMuIndep[i][x][jQ][jPt][n] *                     as * pdflc * BinSize[i];
			  XSectionMuVar_LO[i]	+=  BlockB_LO->SigmaTildeMuFDep [i][x][jQ][jPt][n] * std::log(muf2/Q2) * as * pdflc * BinSize[i];
			  XSectionMuVar_LO[i]	+=  BlockB_LO->SigmaTildeMuRDep [i][x][jQ][jPt][n] * std::log(mur2/Q2) * as * pdflc * BinSize[i];
			}
			// NLO Block
			for(int n=0;n<BlockB_NLO->NSubproc;n++){ 
			  double as	= BlockB_NLO->AlphasTwoPi[i][jQ][jPt];
			  double pdflc	= BlockB_NLO->PdfLcMuVar[i][x][jQ][jPt][n];
			  XSectionMuVar[i]	+=  BlockB_NLO->SigmaTildeMuIndep[i][x][jQ][jPt][n] *                     as * pdflc * BinSize[i];
			  XSectionMuVar[i]	+=  BlockB_NLO->SigmaTildeMuFDep [i][x][jQ][jPt][n] * std::log(muf2/Q2) * as * pdflc * BinSize[i];
			  XSectionMuVar[i]	+=  BlockB_NLO->SigmaTildeMuRDep [i][x][jQ][jPt][n] * std::log(mur2/Q2) * as * pdflc * BinSize[i];
			
			}
		      }
		      
		      //  -----  pp, ppbar ---- //
		      if ( BlockB_LO->NPDFDim == 1 ) {
			// LO Block
			for(int n=0;n<BlockB_LO->NSubproc;n++){ 
			  double as	= BlockB_LO->AlphasTwoPi[i][jQ][jPt];
			  double pdflc	= BlockB_LO->PdfLcMuVar[i][x][jQ][jPt][n];
			  XSectionMuVar[i]	+=  BlockB_LO->SigmaTildeMuIndep[i][x][jQ][jPt][n] *                  as * pdflc * BinSize[i];
			  XSectionMuVar[i]	+=  BlockB_LO->SigmaTildeMuFDep [i][x][jQ][jPt][n] * std::log(muf2) * as * pdflc * BinSize[i];
			  XSectionMuVar[i]	+=  BlockB_LO->SigmaTildeMuRDep [i][x][jQ][jPt][n] * std::log(mur2) * as * pdflc * BinSize[i];

			  XSectionMuVar_LO[i]	+=  BlockB_LO->SigmaTildeMuIndep[i][x][jQ][jPt][n] *                  as * pdflc * BinSize[i];
			  XSectionMuVar_LO[i]	+=  BlockB_LO->SigmaTildeMuFDep [i][x][jQ][jPt][n] * std::log(muf2) * as * pdflc * BinSize[i];
			  XSectionMuVar_LO[i]	+=  BlockB_LO->SigmaTildeMuRDep [i][x][jQ][jPt][n] * std::log(mur2) * as * pdflc * BinSize[i];
			}
			// NLO Block
			for(int n=0;n<BlockB_NLO->NSubproc;n++){ 
			  double as	= BlockB_NLO->AlphasTwoPi[i][jQ][jPt];
			  double pdflc	= BlockB_NLO->PdfLcMuVar[i][x][jQ][jPt][n];
			  XSectionMuVar[i]	+=  BlockB_NLO->SigmaTildeMuIndep[i][x][jQ][jPt][n] *                  as * pdflc * BinSize[i];
			  XSectionMuVar[i]	+=  BlockB_NLO->SigmaTildeMuFDep [i][x][jQ][jPt][n] * std::log(muf2) * as * pdflc * BinSize[i];
			  XSectionMuVar[i]	+=  BlockB_NLO->SigmaTildeMuRDep [i][x][jQ][jPt][n] * std::log(mur2) * as * pdflc * BinSize[i];
			
			}
		      }
		    }
	}
      }
    }
    
  }
}

//______________________________________________________________________________


void FastNLOReader::SetAlphasMz( double AlphasMz ){
  
  if ( AlphasMz != fAlphasMz ){
    fAlphasMz	= AlphasMz;		// new alpha_s value
    FillAlphasCache();
  }
  else {
    // nothing todo!
  }
  
}



//______________________________________________________________________________


void FastNLOReader::FillAlphasCache(){
  
  FillAlphasCacheInBlockB( BlockB_LO  );
  FillAlphasCacheInBlockB( BlockB_NLO );

}



//______________________________________________________________________________


double FastNLOReader::GetAlphas( double Q ){
  //switch (AlphasEvolution )
  if ( fAlphasEvolution == kGRV )			return GetAlphasGRV	( Q , fAlphasMz );
  else if ( fAlphasEvolution == kNLOJET )		return GetAlphasNLOJET	( Q , fAlphasMz );
  else if ( fAlphasEvolution == kCTEQpdf )		return GetAlphasCTEQpdf	( Q , fAlphasMz );
  else if ( fAlphasEvolution == kFastNLO )		return GetAlphasFastNLO	( Q , fAlphasMz );
  else if ( fAlphasEvolution == kLHAPDFInternal )	return GetAlphasLHAPDF	( Q );
  else return 0;
}



//______________________________________________________________________________


void FastNLOReader::FillAlphasCacheInBlockB( FastNLOBlockB* B ){
  

  for(int i=0;i<NObsBin;i++){
    for(int j=0;j<B->GetTotalScalenodes();j++){
      int scalenode1 = j;
      int scalenode2 = j;
      if (B->NScaleDim>1){          
	scalenode1 = j / B->Nscalenode[1];
	     scalenode2 = j % B->Nscalenode[1];
      }
      
      double mur	= B->ScaleNode[i][0][fScalevar][scalenode1];
      double as	= GetAlphas(mur);
      
      double alphastwopi = pow( as/TWOPI , B->Npow );
      B->AlphasTwoPi_v20[i][j] = alphastwopi;
    }
  
    if ( B->NScaleDep == 2 ){
	//printf("fnloBlockB::CalcXsection Xsection2Scales ...");
	for(int jNodeR=0;jNodeR<B->NscalenodeScale1;jNodeR++){
	  double mur		= B->Scale1Node[i][fScalevar][jNodeR];
	  double as		= GetAlphas(mur);
	  double alphastwopi	= pow( as/TWOPI, B->Npow );
	  B->AlphasTwoPi2Scales[i][jNodeR]	= alphastwopi;
	}
    }
    
    if ( B->NScaleDep == 3 ){
	for(int jQ=0;jQ<B->NscalenodeScaleQ;jQ++){
	  for(int jPt=0;jPt<B->NscalenodeScalePt;jPt++){
	    double Q2   = B->ScaleNodeQ[i][jQ]*B->ScaleNodeQ[i][jQ];
	    double Pt   = B->ScaleNodePt[i][jPt];
	    
	    //
	    //   TODO
	    //
	    // todo! here the formula for mu_f and mu_r calculation has to be implemented
	    // grepme choose the scales
	    //      		    double mur2 = (Q2 + Pt*Pt)/2;
	    //     		    double muf2 = mur2;
	    
	    double muf2 = Q2;
	    double as		= GetAlphas(sqrt(muf2));
	    double alphastwopi	= pow( as/TWOPI, B->Npow );
	    B->AlphasTwoPi[i][jQ][jPt] = alphastwopi;
	  }
	}
    }
  }
}


//______________________________________________________________________________


double FastNLOReader::GetAlphasLHAPDF(double Q){
   // the alpha_s evolution as used in LHAPDF
  return LHAPDF::alphasPDF(Q); 
}


//______________________________________________________________________________


double FastNLOReader::GetAlphasNLOJET(double Q, double alphasMZ){

  // this is the evolution, with the one you can make the very very small fnlo error
  // use:
  //   alphasmz = 0.1179 in fnloreader 
  //   double b0  = 1.2202;
  //   double b1  = 0.4897;
  //   double Mz	= 91.70;
  //
  // as evolution motviated by lhpdf.c
  //

  // the original parameters from the cteq6 pdf
  //   // Alpha QCD //
  //   1, 1, 0, -1, 0.1179, 91.70, 1.3, 4.5, 180.0,

  // #define b0 1.2202
  // #define b1 0.4897
  // #define b2 0.1913
  //   double BETA0 =  (11. - 2./3.*NF); // The beta coefficients of the QCD beta function
  //   double BETA1 =  (51. - 19./3.*NF);

  double b0  = 1.2202;
  double b1  = 0.4897;

  //double Mz	= 91.187;
  double Mz	= 91.70;
  //   //   double Mz	= MZ;
  double L = log(Q/Mz);
  //As =  alphasMz; true!
  
   //   L = b0*L;
  L = (b0 + alphasMZ*b1)*L;

   double alphas = alphasMZ/(1.0 + alphasMZ*L);
   //cout << "alphas: " << alphas << "\tQ: " << Q << "\talphasMZ: " << alphasMZ << endl;
   return alphas;
}


//______________________________________________________________________________


double FastNLOReader::GetAlphasGRV(double MU, double ALPSMZ){
  // *******************************************************************
  //    the version from v1.4 as used in fnh2003-desy07073 HERA-I HighQ2 jets
  //    aus: fn-alphas-demo.f 
  // *******************************************************************
  // * M. Wobisch  25/05/99
  // *
  // * calculation of alpha_s in the MSbar scheme
  // * for given alpha_s(Mz)
  // *
  // * using exact / iterative solution of 2-loop formula
  // *
  // * as GRV hep-ph/9806404
  // *
  // * This c++ translation:
  // * D. Britzger  03/03/11
  // *
  // *******************************************************************

  // c - initialize pi and beta functions
  const int NF	= 5;
  const double PI4= TWOPI*2;
  const double B0 =  (11. - 2./3.*NF); // The beta coefficients of the QCD beta function
  const double B1 =  (102. - 38./3.*NF);
  const double B10 = B1 / B0 / B0;
  const double ZMASS2 = pow ( 91.187 , 2 );

  // c - exact formula to extract Lambda from alpha_s(Mz)
  double Q2 = pow(MU,2);
  double LAM2 = ZMASS2 * exp( -1.*PI4/B0/ALPSMZ +  B10 * log( PI4/B0/ALPSMZ + B10) );

  // c - extract approx. alpha_s(mu) value 
  double LQ2 = log( Q2 / LAM2 ) ;
  double ASAPPROX = PI4/B0/LQ2 * (1. - B10*log(LQ2)/LQ2);
  double ALPHAS = ASAPPROX;
    
  // c - exact 2loop value by Newton procedure
  for ( int I = 1 ; I <=6 ; I++ ){
    double  F  = LQ2 - PI4/B0/ALPHAS + B10*log(PI4/B0/ALPHAS + B10);
    double FP = -1.*PI4/B0/(ALPHAS*1.01) + B10 * log(PI4/B0/(ALPHAS*1.01) + B10);
    double FM = -1.*PI4/B0/(ALPHAS*0.99) + B10 * log(PI4/B0/(ALPHAS*0.99) + B10);
    ALPHAS = ALPHAS - F/(FP-FM)*0.02*ALPHAS;
    // c      WRITE(*,*) ' LAMDA/a_s_approx/a_s = ',sqrt(lam2),ASAPPROX,ALPHAS
    // c        alpsmz = alpsmz + I*0.001
    // c       WRITE(*,*) ' alpsmz =', real(alpsmz)
  }
  
  //printf("FNLO v1.4 as.  Q=%7.4f  alphasMZ_0=%7.4f, alphas(Q)=%7.4f\n",Q,alphasMZ,ALPHAS);

  // c - that's it!
  return ALPHAS;

  //       DOUBLE PRECISION FUNCTION ALPS_IT(MU,ALPSMZ,NLOOP)
  // *
  // *  1st version - only for 2-loop - as GRV code
  // *
  //       IMPLICIT NONE
  //       DOUBLE PRECISION MU, ALPSMZ, ALPS4_IT
  //       DOUBLE PRECISION B0, B1, B10 , PI4, F, FP,FM,
  //      +     ONED, TWOD, ZMASS, ZMASS2, ALPHAS, ASAPPROX, Q2, LAM2, LQ2
  //       INTEGER  NLOOP, NF, IFIRST, I, j
  //       PARAMETER (ZMASS = 91.187)        ! PDG data book '98
  // c       PARAMETER (ZMASS = 100)
  //       DATA IFIRST/0/, ONED/1.D0/, TWOD/2.D0/
  //       SAVE IFIRST, NF, ONED, TWOD, PI4, B0, B1, B10, ZMASS2

  // c - initialize pi and beta functions
  //       IF (IFIRST.eq.0) THEN
  //          IFIRST = 1
  // c         WRITE(*,*) '  *   ALPS_IT:  exact 2-loop result for alpha_s'
  //          NF = 5
  //          PI4 = 4D0 * 4D0 * ATAN(1D0)
  //          B0  = 11D0 - 2D0/3D0 * DBLE(NF)
  //          B1  = 102D0 - 38D0 / 3D0 * DBLE(NF)
  //          B10 = B1 / B0 / B0
  //          ZMASS2 = ZMASS**2
  //          IF (NLOOP .ne. 2) WRITE(*,*) 'ALPS_IT:  only for 2-loop!!'
  //       ENDIF

  // c - exact formula to extract Lambda from alpha_s(Mz)
  //       Q2 = MU**2
  //       LAM2 = ZMASS2 * EXP( -PI4/B0/(ALPSMZ) + 
  //      +     B10 * DLOG( PI4/B0/ALPSMZ + B10) )

  // c - extract approx. alpha_s(mu) value 
  //       LQ2 = DLOG( Q2 / LAM2 ) 
  //       ASAPPROX = PI4/B0/LQ2 * (1D0 - B10*DLOG(LQ2)/LQ2)
  //       ALPHAS = ASAPPROX

  // c - exact 2loop value by Newton procedure
  //       DO I=1,6
  //          F  = LQ2 - PI4/B0/ALPHAS + B10*DLOG(PI4/B0/ALPHAS + B10)
  //          FP = - PI4/B0/(ALPHAS*1.01D0) + 
  //      +        B10 * DLOG(PI4/B0/(ALPHAS*1.01D0) + B10)
  //          FM = - PI4/B0/(ALPHAS*0.99D0) + 
  //      +        B10 * DLOG(PI4/B0/(ALPHAS*0.99D0) + B10)
  //          ALPHAS = ALPHAS - F/(FP-FM)*0.02D0*ALPHAS 
  // c      WRITE(*,*) ' LAMDA/a_s_approx/a_s = ',sqrt(lam2),ASAPPROX,ALPHAS
  // c        alpsmz = alpsmz + I*0.001
  // c       WRITE(*,*) ' alpsmz =', real(alpsmz)
  //       ENDDO

  // c - that's it!
  //       ALPS_IT = ALPHAS
  //       RETURN 
  //       END
  // }

}


//______________________________________________________________________________


double FastNLOReader::GetAlphasCTEQpdf(double Q, double alphasMZ){

  //as evolution motivated by cteq-pdf-1.0.4 as used in nlojet 4.1.3 without crack

  double as_twopi = alphasMZ/TWOPI;
  double Mz	= 91.187;
  //double Mz	= 91.70;
  //   double Mz	= MZ;

  int ord=2;
  int nf=5;

  const int NF	= 5;

  double b0 =  (11. - 2./3.*NF); // The beta coefficients of the QCD beta function
  double b1 =  (51. - 19./3.*NF);

  //double astolmd(unsigned int ord, unsigned int nf, double as, double q)
  // ----------------------------------
  double t8 = 1.0/(b0*as_twopi);
  //       /*  it is a leading order evolution  */
  //       if(ord <= 1) return q*exp(-t8);
  /*  at NLO or higer order level returns with the NLO alpha_s  */
  double as0, as1, ot, lt, br = (51.0 - 19.0/3.0*nf)/(b0*b0);
  do {
    lt = log(2.0*t8)/t8;
    ot = t8;

    as0 = (1.0 - br*lt)/(b0*t8);
    as1 = (-1.0 - br*(1.0/t8-2.0*lt))/(b0*t8*t8);
    t8 += (as_twopi - as0)/as1;
  } while(fabs(ot-t8)/ot > 1e-5);
  double lmd = Mz*exp(-t8);
  // ----------------------------------

  double t = log(Q/lmd);
  double asMz = 1.0/(b0*t);

  return asMz*(1.0-b1/b0*asMz*log(2.0*t)) *TWOPI;

}



//______________________________________________________________________________


double FastNLOReader::GetAlphasFastNLO(double Q, double alphasMZ){
  //
  //  the original FNLOv2.0 implementation
  //
  //   cout << "using old alphas evolution." << endl;
  const int NF	= 5;
  double BETA0 =  (11. - 2./3.*NF); // The beta coefficients of the QCD beta function
  double BETA1 =  (51. - 19./3.*NF);

  //    // This is from NLOJET++, alpha.cc
  double Mz	= 91.187;
  double res	= alphasMZ;
  double b0	= BETA0/TWOPI;
  double w = 1.0 + b0*alphasMZ*log(Q/Mz);
  res /= w;
  double b1 = BETA1/TWOPISQR;
  res *= 1.0 - alphasMZ*b1/b0*log(w)/w;
  return res;

}




//______________________________________________________________________________


void FastNLOReader::FillPDFCache(){
  
  if ( fLHAPDFfilename == "" || fLHAPDFpath ==""){
    printf("ERROR. You must specify a LHAPDF filename first.\n"); exit(1);
  }


  if ( fPDFInterface == kLHAPDF ){
    InitLHAPDF();
    FillBlockBPDFLCsWithLHAPDF(BlockB_LO);
    FillBlockBPDFLCsWithLHAPDF(BlockB_NLO);
  }
  else if ( fPDFInterface == kH1FITTER ){
    FillBlockBPDFLCsWithH1Fitter(BlockB_LO);
    FillBlockBPDFLCsWithH1Fitter(BlockB_NLO);
  }
  
}



//______________________________________________________________________________


void FastNLOReader::InitLHAPDF(){
  
  
  if ( fLHAPDFfilename == "" || fLHAPDFpath ==""){
    printf("FastNLOReader::FillPDFCacheLHAPDF(). ERROR. You must specify a LHAPDF filename first.\n"); exit(1);
  }

  // Do I need to specify the lhapdf path here?? 
  // isn't this done through the linkin??

  string LHAPDFfile = fLHAPDFpath+"/"+fLHAPDFfilename;
  
  // ---- check if file exists ----- //
  FILE* fp = fopen(LHAPDFfile.c_str(), "r");
  if (fp) {
    fclose(fp);
  } else {
    printf("Error. File does not exists. Was looking for: %s. Exiting.\n",LHAPDFfile.c_str());
    exit(1);
  } 

  LHAPDF::setVerbosity(LHAPDF::SILENT);
  LHAPDF::initPDFSetByName(fLHAPDFfilename);
  // lastusedpdf=1; // if you use multiple pdfs and want to switch between them...
  fnPDFs = LHAPDF::numberPDF();
  if ( fnPDFs < fiPDFSet ){
    cout << "Error. There are only " << fnPDFs << " pdf sets within this LHAPDF file. You were looking for set number " << fiPDFSet << endl;
  }

  LHAPDF::initPDF(fiPDFSet);

}



//______________________________________________________________________________


void FastNLOReader::FillBlockBPDFLCsWithH1Fitter( FastNLOBlockB* B ){
  // 
  // here is the inteface of h1fitter with FastNLO.
  // 
  // There are those "BlockB's" which hold the interpolated table, as well
  // as pre-sized tables for the corresponding pdf-linear combinations.
  //
  // 


  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //
  //         THIS CODE IS ONLY A STUMP  !!!!!
  // 
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   // linear: DIS-case
   if(B->NPDFDim == 0){
     vector<double> xfx; // PDFs of all 13 partons
     xfx.resize(13);
     
     for(int i=0;i<NObsBin;i++){				// loop over all observation bins.
       int nxmax = B->GetNxmax(i);		
       for(int j=0;j<B->Nscalenode[0];j++){			// loop over all 'scale'-nodes. Those are nodes of Mu_f in each ObsBin.
	 for(int k=0;k<nxmax;k++){				// loop over all 'x'-nodes where you want to determine your pdf.
		      
	   double xp	= B->XNode1[i][k];			// calculate x
	   double muf	= B->ScaleNode[i][0][fScalevar][j];	// calculate muf

	   // Krzys: you have to change this line!
	   xfx = LHAPDF::xfx(xp,muf);				// get 13 parton distributions here!

	   vector < double > buffer = CalcPDFLinearComb(xfx,xfx,B->IPDFdef1, B->IPDFdef2, B->NSubproc ); //calculate linear combinations as used in fastNLO
	   for(int l=0;l<B->NSubproc;l++){ 
	     B->PdfLc[i][j][k][l] = buffer[l];			// fill the pdf cache within BlockB here!
	   }
	 }
       }
     }
		
     // ---- we also fill the pdfLc for 2Scale interpolation ---- //
     if ( B->NScaleDep == 2 ){
		xfx.clear();
		xfx.resize(13);
		
		for(int i=0;i<NObsBin;i++){
		  int nxmax = B->GetNxmax(i);
		  for(int j=0;j<B->NscalenodeScale2;j++){
		    for(int k=0;k<nxmax;k++){ 
		      double xp		= B->XNode1[i][k];
		      double muf	= B->Scale2Node[i][fScalevar][j];

		      // Krzys: you have to change this line!
		      xfx = LHAPDF::xfx(xp,muf); // LHAPDF::xfx_p_(x,muf,0,0)
		      
		      vector < double > buffer  = CalcPDFLinearComb(xfx,xfx,B->IPDFdef1, B->IPDFdef2, B->NSubproc );
		      for(int l=0;l<B->NSubproc;l++){ 
			B->PdfLc2Scales[i][j][k][l] = buffer[l];
		      }
		    }
		  }
		}
     }
  
   if ( B->NScaleDep == 3 ){
		for(int i=0;i<NObsBin;i++){
		  int nxmax = B->GetNxmax(i);
		  for(int jQ=0;jQ<B->NscalenodeScaleQ;jQ++){
		    for(int jPt=0;jPt<B->NscalenodeScalePt;jPt++){
		      
		      //
		      //   TODO
		      //
		      // todo! here the formula for mu_f calculation has to be implemented
		      // grepme choose the scales
		      //double muf2 = (ScaleNodeQ[i][jQ]*ScaleNodeQ[i][jQ] + ScaleNodePt[i][jPt]*ScaleNodePt[i][jPt] )/ 2.;
		      double muf2 = B->ScaleNodeQ[i][jQ]*B->ScaleNodeQ[i][jQ];
		      
		      for(int x=0;x<nxmax;x++){ 
			double xp	= B->XNode1[i][x];
			
			// Krzys: you have to change this line!
			xfx = LHAPDF::xfx(xp,sqrt(muf2)); // LHAPDF::xfx_p_(x,muf,0,0)
			
			vector < double > buffer  = CalcPDFLinearComb(xfx,xfx,B->IPDFdef1, B->IPDFdef2, B->NSubproc );
			for(int l=0;l<B->NSubproc;l++){ 
			  B->PdfLcMuVar[i][x][jQ][jPt][l] = buffer[l];
			}
		      }
		    }
		  }
		}
		
     }
   }
}



//______________________________________________________________________________


void FastNLOReader::FillBlockBPDFLCsWithLHAPDF( FastNLOBlockB* B ){
  
  //
  // this method fills the already 'resized'
  // pdflc-vectors inside FastNLOBlockB class.
  //

   if (B->NScaleDim>1){
     printf("FastNLOReader::FillBlockBPDFLCsWithLHAPDF. WOW! NScaleDim>1! This is usually not the case!\n");
     //scaleindex2 = 1; // If we use multiple scales, then mu_f is by convention the second scale -> index=1 
     //fScalevar2 = fScalevar % NfScalevar[1]; 
   }


   // linear: DIS-case
   if(B->NPDFDim == 0){
     vector<double> xfx; // PDFs of all partons
     xfx.resize(13);
     
     for(int i=0;i<NObsBin;i++){
       int nxmax = B->GetNxmax(i);
       for(int j=0;j<B->Nscalenode[0];j++){
	 for(int k=0;k<nxmax;k++){ 
		      
	   double xp	= B->XNode1[i][k];
	   double muf	= B->ScaleNode[i][0][fScalevar][j];
	   xfx = LHAPDF::xfx(xp,muf); // LHAPDF::xfx_p_(x,muf,0,0)
	   vector < double > buffer = CalcPDFLinearComb(xfx,xfx,B->IPDFdef1, B->IPDFdef2, B->NSubproc ); //calculate linear combinations
	   for(int l=0;l<B->NSubproc;l++){ 
	     B->PdfLc[i][j][k][l] = buffer[l];
	   }
	 }
       }
     }
		
     // ---- we also fill the pdfLc for 2Scale interpolation ---- //
     if ( B->NScaleDep == 2 ){
		
		for(int i=0;i<NObsBin;i++){
		  int nxmax = B->GetNxmax(i);
		  for(int j=0;j<B->NscalenodeScale2;j++){
		    for(int k=0;k<nxmax;k++){ 
		      double xp		= B->XNode1[i][k];
		      double muf	= B->Scale2Node[i][fScalevar][j];
		      xfx = LHAPDF::xfx(xp,muf); // LHAPDF::xfx_p_(x,muf,0,0)
		      vector < double > buffer  = CalcPDFLinearComb(xfx,xfx,B->IPDFdef1, B->IPDFdef2, B->NSubproc );
		      for(int l=0;l<B->NSubproc;l++){ 
			B->PdfLc2Scales[i][j][k][l] = buffer[l];
		      }
		    }
		  }
		}
     }

     if ( B->NScaleDep == 3 ){
		xfx.clear();
		xfx.resize(13);
		for(int i=0;i<NObsBin;i++){
		  int nxmax = B->GetNxmax(i);
		  for(int jQ=0;jQ<B->NscalenodeScaleQ;jQ++){
		    for(int jPt=0;jPt<B->NscalenodeScalePt;jPt++){
		      
		      //
		      //   TODO
		      //
		      // todo! here the formula for mu_f calculation has to be implemented
		      // grepme choose the scales
		      //double muf2 = (ScaleNodeQ[i][jQ]*ScaleNodeQ[i][jQ] + ScaleNodePt[i][jPt]*ScaleNodePt[i][jPt] )/ 2.;
		      double muf2 = B->ScaleNodeQ[i][jQ]*B->ScaleNodeQ[i][jQ];
		      
		      for(int x=0;x<nxmax;x++){ 
			double xp	= B->XNode1[i][x];
			xfx = LHAPDF::xfx(xp,sqrt(muf2)); // LHAPDF::xfx_p_(x,muf,0,0)
			vector < double > buffer  = CalcPDFLinearComb(xfx,xfx,B->IPDFdef1, B->IPDFdef2, B->NSubproc );
			for(int l=0;l<B->NSubproc;l++){ 
			  B->PdfLcMuVar[i][x][jQ][jPt][l] = buffer[l];
			}
		      }
		    }
		  }
		}
		
     }
   }

   // half matrix (pp, ppbar)
   else if( B->NPDFDim == 1){
     int scaleindex2 = 0;

     if ( B->NScaleDep != 3 ){
       printf("ERROR in fnloBlockB::FillPDFCache. Reading of pp and ppbar tables is only implemented for MuVar tables.\n");
       exit(1);
     }

     
     if ( B->NScaleDep == 3 ){
       vector < vector < double > > xfx; // PDFs of all partons
       for(int i=0;i<NObsBin;i++){
         int nxmax = B->GetNxmax(i);  // total entries in half matrix
         int nxbins1 = B->Nxtot1[i]; // number of columns in half matrix

         xfx.resize(nxbins1);

         //for(int j=0;jB-><Nscalenode[scaleindex2];j++){

	 for(int jQ=0;jQ<B->NscalenodeScaleQ;jQ++){
	   for(int jPt=0;jPt<B->NscalenodeScalePt;jPt++){
		      
	     // determine all pdfs of hadron1
	     for(int k=0;k<nxbins1;k++){ 
	       xfx[k].resize(13);

	       double muf2 = B->ScaleNodeQ[i][jQ]*B->ScaleNodeQ[i][jQ];
	       double xp	= B->XNode1[i][k];
	       xfx[k] = LHAPDF::xfx(xp,sqrt(muf2)); // LHAPDF::xfx_p_(x,muf,0,0)
	     
	     }
	     int x1bin = 0;
	     int x2bin = 0;
	     for(int k=0;k<nxmax;k++){ 
	       
	       vector < double > buffer;
	       if ( B->NPDF != 2 ) exit(1);
	       if ( B->NPDFPDG[0] == B->NPDFPDG[1] )		buffer = CalcPDFLinearCombPPMuVar(xfx[x1bin],xfx[x2bin] );
	       else if ( B->NPDFPDG[0] == -B->NPDFPDG[1] )	buffer = CalcPDFLinearCombPPbarMuVar(xfx[x1bin],xfx[x2bin] );
	       else exit(1); // unknown parton configuration.

	       //vector < double > buffer = CalcPDFLinearComb(xfx[x1bin],xfx[x2bin], B->IPDFdef1, B->IPDFdef2, B->NSubproc); //calculate linear combinations
		 
	       for(int l=0;l<7;l++){  // in MuVar tables, we always have 7 subprocess for pp or ppbar processes, as it is the case in nlojet++
		 B->PdfLcMuVar[i][k][jQ][jPt][l] = buffer[l];
	       }
	       x1bin++;
	       if(x1bin>x2bin){
		 x1bin = 0;
		 x2bin++;
	       }
	     }
	   }
         }
       }
     }
   }



}



//______________________________________________________________________________


vector<double> FastNLOReader::CalcPDFLinearComb(vector<double> pdfx1, vector<double> pdfx2, int IPDFdef1, int IPDFdef2, int NSubproc){
  
     
   switch(IPDFdef1){
   case 2: // ep , DIS and gP
      switch(IPDFdef2){
      case 1: // DIS: determine gluon,sigma,delta,
      case 2: // direct gammaP: gluon,sigma,delta


	return CalcPDFLinearCombDIS( pdfx1 , NSubproc );
	
	break;
      default: printf("fnloBlockB::CalcPDFLinearComb : Ipdfdef1=2, Ipdfdef2= %d not supported. Exit.\n",IPDFdef2); exit(1);
      }
      break;
   case 3:
      switch(IPDFdef2){
      case 1: // ppbar: gg   qg   gq   qr   qq   qqb   qrb 

	printf("FastNLOReader::CalcPDFLinearComb. Please call CalcPDFLinearCombPPMuVar or CalcPDFLinearCombPPbarMuVar directly.\n");
	return CalcPDFLinearCombPPbarMuVar( pdfx1 , pdfx2 );

	/*
	// ---------------------------------------------------
	// remember from nlojet++ proc-hhc/weight.cc
	//   const char *weight_label_hhc[7] = {"gg", "qg", "gq", "qr", "qq", "qqb", "qrb"};
	// ---------------------------------------------------
	// remember from nlojet++ proc-hhc/process.cc
	//    retval[0] = A0*B0;
	//    retval[1] = (A + Ab)*B0;
	//    retval[2] = A0*(B + Bb);
	//    retval[3] = A*B + Ab*Bb - D;
	//    retval[4] = D;
	//    retval[5] = Db;
	//    retval[6] = A*Bb +Ab*B - Db;
	// ---------------------------------------------------
 

	//          double B0,B,Bb;
	//          double A0,A,Ab;
	//          double D,Db;

	//          A0 = pdfx2[6];
	//          B0 = pdfx1[6];
         
	//          A = Ab = 0.;
	//          B = Bb = 0.;
	//          D = Db = 0.;
	//           for(int l=0;l<6;l++){
	//              A  += pdfx2[l+7];
	//              Ab += pdfx2[5-l];
	// // pp
	// //              B  += pdfx1[l+7];
	// //              Bb += pdfx1[5-l];
	// //              D  += pdfx1[l+7] * pdfx2[l+7] + pdfx1[5-l] * pdfx2[5-l];
	// //              Db += pdfx1[l+7] * pdfx2[5-l] + pdfx1[5-l] * pdfx2[l+7];
	// // ppbar
	//              Bb  += pdfx1[l+7];
	//              B += pdfx1[5-l];
	//              Db  += pdfx1[l+7] * pdfx2[l+7] + pdfx1[5-l] * pdfx2[5-l];
	//              D += pdfx1[l+7] * pdfx2[5-l] + pdfx1[5-l] * pdfx2[l+7];
	//           }         

	//          pdflc[0] = A0*B0; // gluon gluon
	//          pdflc[1] = (A + Ab)*B0; // quark gluon
	//          pdflc[2] = A0*(B + Bb);
	//          pdflc[3] = A*B + Ab*Bb - D;
	//          pdflc[4] = D;
	//          pdflc[5] = Db;
	//          pdflc[6] = A*Bb +Ab*B - Db;        

	 
	// 	working code from example.f
	// * input:
	// *    ireact             flag for reaction (1:DIS, 2:pp-jets, 3:ppbar-jets)
	// *    i                  x-index of first hadron
	// *    j                  x-index for second hadron (if two-hadron process)
	// *    XPDF1(nxmax,-6:6)  PDF array for all x-bins
	// *    XPDF2(nxmax,-6:6)  PDF array for all x-bins
	// *
	// * output:
	// *    H(10)              PDF linear combinations
	// *
	//       Double Precision XPDF1(MxNxTot,-6:6),XPDF2(MxNxTot,-6:6), H(10),
	//      +     G1, G2,              ! gluon densities from both hadrons
	//      +     SumQ1, SumQ2,        ! sum of quark densities
	//      +     SumQB1, SumQB2,      ! sum of anti-quark densities
	//      +     Q1(6),Q2(6), QB1(6),QB2(6), ! arrays of 6 (anti-)quark densities
	//      +     S,A                  ! products S,A

	// 	c --- hadron-hadron: jets
	// 	Elseif (icf1.eq.3.and.icf2.eq.1.and.(icf3.ge.1.and.icf3.le.2))Then 
	 
	//  	 double G1  = 0;
	//  	 double G2  = 0;
	//  	 double SumQ1  = 0;
	//  	 double SumQB1 = 0;
	//  	 double SumQ2  = 0;
	//  	 double SumQB2 = 0;
	//  	 double Q1[6] = {0};
	//  	 double QB1[6] = {0};
	//  	 double Q2[6] = {0};
	//  	 double QB2[6] = {0};
	//  	 double S = 0;
	//  	 double A = 0;

	double G1 ;
	double G2 ;
	double SumQ1;
	double SumQB1;
	double SumQ2 ;
	double SumQB2;
	double Q1[6] ;
	double QB1[6];
	double Q2[6] ;
	double QB2[6];
	double S ;
	double A ;

	G1  = 0;
	G2  = 0;
	SumQ1  = 0;
	SumQB1 = 0;
	SumQ2  = 0;
	SumQB2 = 0;
	S = 0;
	A = 0;

	// *** Here ***
	// indices for pdfx1 and pdfx2:
	// 0..5 = tbar, ..., ubar, dbar;
	// 6 = g;
	// 7..12 = d, u, ..., t
	// 
	// *** fortran reader ***
	// -6 - -1
	// 0 gluon
	// 1 - 6
	//
	//  =>  k_fortran -> k+6 , when accessing the pdf

	for ( int k = 0 ; k < 6 ; k++ ){
	//	   Do k=1,6
	Q1[k]  = pdfx1[k+6];//  XPDF1(i,k)  ! read 1st PDF at x1
	QB1[k] = pdfx1[(-k)+6];//XPDF1(i,-k)
	SumQ1  += Q1[k];
	SumQB1 += QB1[k];
	Q2[k]  = pdfx2[k+6];//XPDF2(j,k)  ! read 2nd PDF at x2
	QB2[k] = pdfx2[(-k)+6];//XPDF2(j,-k)
	SumQ2  +=  Q2[k];
	SumQB2 += QB2[k];
	//Enddo
	}
	G1     =  pdfx1[6];//XPDF1(i,0)
	G2     =  pdfx2[6];//XPDF1(j,0)
	//c   - compute S,A
	for ( int k = 0 ; k < 6 ; k++ ){
	//Do k=1,6
	S += (Q1[k]*Q2[k]) + (QB1[k]*QB2[k]); 
	A += (Q1[k]*QB2[k]) + (QB1[k]*Q2[k]); 
	//  Enddo
	}
	//c   - compute seven combinations

	pdflc[0] = G1*G2;
	pdflc[1] = SumQ1*SumQ2 + SumQB1*SumQB2 - S;
	pdflc[2] = S;
	pdflc[3] = A;
	pdflc[4] = SumQ1*SumQB2 + SumQB1*SumQ2 - A;
	pdflc[5] = (SumQ1+SumQB1)*G2;
	pdflc[6] = G1*(SumQ2+SumQB2);
	 
	// if 6 subprocesses: sum pdflc[6]+pdflc[7] in the calling function
	//If (icf3.eq.1) H[6] = H[6]+H[7] ! case: 6 subproc

	break;
	default: printf("fnloBlockB::CalcPDFLinearComb :Ipdfdef1=3, Ipdfdef2= %d not supported. Exit.\n",IPDFdef2); exit(1);
	}
	break;
	case 4:
	switch(IPDFdef2){
	case 1: // resolved gammaP: gg   qg   gq   qr   qq   qqb   qrb 
	double B0,B,Bb;
	double A0,A,Ab;
	double D,Db;

	A0 = pdfx2[6];
	B0 = pdfx1[6];
         
	A = Ab = 0.;
	B = Bb = 0.;
	D = Db = 0.;
	for(int l=0;l<6;l++){
	A  += pdfx2[l+7];
	Ab += pdfx2[5-l];
	B  += pdfx1[l+7];
	Bb += pdfx1[5-l];
	D  += pdfx1[l+7] * pdfx2[l+7] + pdfx1[5-l] * pdfx2[5-l];
	Db += pdfx1[l+7] * pdfx2[5-l] + pdfx1[5-l] * pdfx2[l+7];
	}         

	pdflc[0] = A0*B0; // gluon gluon
	pdflc[1] = (A + Ab)*B0; // quark gluon
	pdflc[2] = A0*(B + Bb);
	pdflc[3] = A*B + Ab*Bb - D;
	pdflc[4] = D;
	pdflc[5] = Db;
	pdflc[6] = A*Bb +Ab*B - Db;
	*/ 
	break;
      default: printf("fnloBlockB::CalcPDFLinearComb ::Ipdfdef1=4, Ipdfdef2= %d not supported. Exit.\n",IPDFdef2); exit(1);
      }
      break;
   default: printf("fnloBlockB::CalcPDFLinearComb :Ipdfdef1= %d not supported. Exit.\n",IPDFdef1); exit(1);
   }

  
}



//______________________________________________________________________________


vector<double> FastNLOReader::CalcPDFLinearCombPPbarMuVar(vector<double> pdfx1, vector<double> pdfx2 ){
  
  // just 'flip' pdfx2 and call PP-LC
  vector < double > pdfx2bar;
  pdfx2bar.resize(pdfx2.size());
  if ( pdfx2bar.size() != 13 ){
    printf("FastNLOReader::CalcPDFLinearCombPPbarMuVar. Assuming that pdf has 13 partons but current pdf has %d.\n",pdfx2bar.size());
  }
  for ( int i = 0 ; i < 13 ; i++ ){
    pdfx2bar[i]	= pdfx2[12-i];
  }
  
  return CalcPDFLinearCombPPMuVar(pdfx1,pdfx2bar);
    
}



//______________________________________________________________________________


vector<double> FastNLOReader::CalcPDFLinearCombPPMuVar(vector<double> pdfx1, vector<double> pdfx2 ){
 
  vector < double > retval;
  retval.resize(7);
  
  // ---------------------------------------------------
  // remember from nlojet++ proc-hhc/weight.cc
  //   const char *weight_label_hhc[7] = {"gg", "qg", "gq", "qr", "qq", "qqb", "qrb"};
  // ---------------------------------------------------
  // remember from nlojet++ proc-hhc/process.cc
  //    retval[0] = A0*B0;
  //    retval[1] = (A + Ab)*B0;
  //    retval[2] = A0*(B + Bb);
  //    retval[3] = A*B + Ab*Bb - D;
  //    retval[4] = D;
  //    retval[5] = Db;
  //    retval[6] = A*Bb +Ab*B - Db;
  // ---------------------------------------------------
  unsigned int nu	= 2;
  unsigned int nd	= 3;

  int ia, iq;
  //weight_hhc retval;
  //   static double __f1[13], __f2[13];
  static const int iu[3] = {2,4,6}, id[3] = {1,3,5};
  
  //   //----- calculat the pdfs -----
  //   double *f1 = __f1+6, *f2 = __f2+6;
  
  //   this -> hadronA(x1, mf2, nu, nd, f1);
  //   this -> hadronB(x2, mf2, nu, nd, f2);
  
  //----- gluon pdfs -----
  double A0 = pdfx1[0+6];
  double B0 = pdfx2[0+6];
  
  
  //---- up type quarks -----
  double q1, q2, a1, a2;
  double A = 0.0, B = 0.0, Ab = 0.0, Bb = 0.0, D = 0.0, Db = 0.0; 
  
  for(unsigned int u = 0; u < nu && u < 3; u++) {
    ia = -(iq = iu[u]);
    q1 = pdfx1[iq+6]; q2 = pdfx2[iq+6];
    a1 = pdfx1[ia+6]; a2 = pdfx2[ia+6];
      
    A += q1; Ab += a1; B += q2; Bb += a2;
    D += q1*q2 + a1*a2; Db += q1*a2 + a1*q2;
  }
    
  //----- down type quarks -----
  for(unsigned int d = 0; d < nd && d < 3; d++) {
    ia = -(iq = id[d]);
    //cout << "ia="<<ia<<"\tiq="<<iq<<"\td="<<d<< endl;
    q1 = pdfx1[iq+6]; q2 = pdfx2[iq+6];
    a1 = pdfx1[ia+6]; a2 = pdfx2[ia+6];
      
    A += q1; Ab += a1; B += q2; Bb += a2;
    D += q1*q2 + a1*a2; Db += q1*a2 + a1*q2;
  }
    
  retval[0] = A0*B0;
  retval[1] = (A + Ab)*B0;
  retval[2] = A0*(B + Bb);
  retval[3] = A*B + Ab*Bb - D;
  retval[4] = D;
  retval[5] = Db;
  retval[6] = A*Bb +Ab*B - Db;

  return retval;
 
}



//______________________________________________________________________________


vector<double> FastNLOReader::CalcPDFLinearCombDIS( vector<double> pdfx1 , int NSubproc){
  
  vector < double > pdflc;
  pdflc.resize(3);

  pdflc[1] = pdfx1[6]; //gluon

  
  for(int l=0;l<13;l++){
    double temp = (l==6 ? 0.0 : pdfx1[l]);
    if (!(l&1)) temp *= 4.;
    pdflc[0] += temp; // delta
  }
  pdflc[0] /= 9.;
  
  if(NSubproc>2){ // only from NLO
    for(int l=0;l<6;l++){
      pdflc[2] += pdfx1[5-l] + pdfx1[l+7]; // sigma
    }
  }

  return pdflc;
    
}



//______________________________________________________________________________


