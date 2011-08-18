// Author: Daniel Britzger
// DESY, 23/07/2011

//  Version 0.4, 
//
//  History:
//    Version 0, initial version

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  FastNLOReader                                                       //
//                                                                      //
//  FastNLOReader is a standalone code for reading                      //
//  FastNLO tables of version 2.0 for DIS processes                     //
//  It is also optimized for an integration into                        //
//  the H1Fitter project.                                               //
//                                                                      //
//  FastNLO is developed by                                             //
//    D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch         //
//    (publication in preparation)                                      //
//    http://projects.hepforge.org/fastnlo                              //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "FastNLOReader.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include "LHAPDF.h"

using namespace std;


//______________________________________________________________________________


extern "C"{
  void fpdfxq_(int *iset, const double *x, const double *q2, double *pdfs,int *ichk);
  void evolution_();
  //int getalf_( double* alfs, double* r2 );
}


//______________________________________________________________________________


FastNLOReader::FastNLOReader(void)
{
   //
   // do not call the standard constructor
   // 
   BlockB_LO		= NULL;
   BlockB_NLO		= NULL;
   BlockB_LO_Ref	= NULL;
   BlockB_NLO_Ref	= NULL;
   cout << "FastNLOReader::FastNLOReader. Please set a filename using SetFilename(<name>)! "<<endl;
}


FastNLOReader::FastNLOReader(string filename)
{
   BlockB_LO		= NULL;
   BlockB_NLO		= NULL;
   BlockB_LO_Ref	= NULL;
   BlockB_NLO_Ref	= NULL;
   
   SetFilename(filename);
}


//______________________________________________________________________________


FastNLOReader::~FastNLOReader(void)
{
   if ( BlockB_LO )		delete BlockB_LO;
   if ( BlockB_NLO )		delete BlockB_NLO;
   if ( BlockB_LO_Ref )		delete BlockB_LO_Ref;
   if ( BlockB_NLO_Ref )	delete BlockB_NLO_Ref;
}


//______________________________________________________________________________



void FastNLOReader::SetFilename(string filename){
   ffilename	= filename;
   Init();
}


//______________________________________________________________________________



void FastNLOReader::Init(){
  printf(" ***************************************************************** \n");
  printf(" *  \n");
  printf(" *  FastNLO Reader - version 0.4\n");
  printf(" *  \n");
  printf(" *  This code is a reader for FastNLO tables, which were\n");
  printf(" *  calculated using nlojet++ v4.1.3 and FastNLO v2.0 (or higher).\n"); 
  printf(" *  \n"); 
  printf(" *  The FastNLO package by\n"); 
  printf(" *    D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch\n");
  printf(" *    (publication in preparation)\n");
  printf(" *  \n");
  printf(" *  for NLOJET++ please cite \n");
  printf(" *    Z. Nagy, Phys. Rev. Lett. 88, 122003 (2002),\n");
  printf(" *    Z. Nagy, Phys. Rev. D68, 094002 (2003)\n");
  printf(" *  \n");
  printf(" ***************************************************************** \n");

  ReadTable();

  SetPDFInterface(FastNLOReader::kLHAPDF);
  SetAlphasEvolution(FastNLOReader::kGRV);
  
  InitScalevariation();

}


//______________________________________________________________________________



void FastNLOReader::InitScalevariation(){
  
  fScaleFacMuR	= 1.;
  fScaleFacMuF	= 1.;
  fScalevar	= 0;

  if ( BlockB_NLO->NScaleDep == 0 ){
    // this is an 'original' v2.0 table
    printf (" *  This table has following %d scale variations for 'theory-error' determination.\n",BlockB_NLO->Nscalevar[0]);
    printf (" *    scalevar #n -> scalefactor\n");
    for ( int i = 0 ; i<BlockB_NLO->Nscalevar[0]; i++ ){
      printf (" *         '%d'    ->    %4.2f .\n", i, BlockB_NLO->ScaleFac[0][i]);
    }
    printf (" *    Setting scale factor to %4.2f and varying mu_f and mu_r simultaneously.\n",BlockB_NLO->ScaleFac[0][0]);
    fScalevar	= 0;
  }

  else if ( BlockB_NLO->NScaleDep == 3 ){
    // this is a MuVar table. You can vary mu_f and mu_r independently by any factor
    // and you can choose the functional form of mu_f and mu_r as functions of
    // scale1 and scale1 (called partly scaleQ2 and scalePt).
    
    if ( BlockB_NLO->NscaleDescript[0] <0 ) {
      printf("Error. No scaledescription available.\n"); // the code will crash soon.
      fMuFFunc	= kScale1;
      fMuRFunc	= kScale1;
      return;
    }

    // ---- DIS ---- //
    if ( BlockB_LO->NPDFDim == 0 ) {
      fMuRFunc	= kQuadraticMean;
      fMuFFunc	= kScale1;
      printf (" *    Setting factorization scale to mu_f^2 = %s^2 .\n", BlockB_NLO->ScaleDescript[0][0].c_str() );
      if ( BlockB_NLO->NscaleDescript[0] == 2 ){
	printf (" *    Setting renormalization scale to mu_r^2 = (%s^2 + %s^2)/2 .\n", BlockB_NLO->ScaleDescript[0][0].c_str() , BlockB_NLO->ScaleDescript[0][1].c_str() );
      }
      else if ( BlockB_NLO->NscaleDescript[0] == 1 &&  BlockB_LO->NscalenodeScalePt > 3 ){
	printf("FastNLOReader::InitScalevariation. Warning. Could not find description for scale variables.\n");
	printf (" *    Setting renormalization scale to mu_r^2 = (scale1^2 + scale2^2)/2 .\n" );
      }
      else if ( BlockB_NLO->NscaleDescript[0] == 1 &&  BlockB_LO->NscalenodeScalePt <= 3  ){
	printf("FastNLOReader::InitScalevariation. Warning. This table has only one scale variable %s stored.\n", BlockB_NLO->ScaleDescript[0][0].c_str() );
	printf (" *    Setting renormalization scale to mu_r^2 = %s^2 .\n", BlockB_NLO->ScaleDescript[0][0].c_str() );
	fMuRFunc	= kScale1;
      }
      else {
	printf("Error. I don't know what to do.\n");
      }
    }
        
    else printf("Error. Unknown process.\n");

  }
  
  else {
    printf("FastNLOReader::InitScalevariation(). ERROR. Could not identify table..\n");
  }
  
  
}


//______________________________________________________________________________



double FastNLOReader::CalcMu( FastNLOReader::EMuX kMuX , double scale1, double scale2, double scalefac ){
  //
  //  Calculate the scales with the defined function and the 
  //  corresponding prefactor.
  //
  
  if ( kMuX == kMuR && fScaleFacMuR != scalefac ) printf("Error. Sth. went wrong with the scales.\n");
  if ( kMuX == kMuF && fScaleFacMuF != scalefac ) printf("Error. Sth. went wrong with the scales.\n");
  
  EScaleFunctionalForm Func;
  if ( kMuX  == FastNLOReader::kMuR )		Func	= fMuRFunc;    // return renormalization scale
  else						Func	= fMuFFunc;    // return factorization scale
  //   else if ( kMuX  == FastNLOReader::kMuF )	Func	= fMuFFunc;    // return factorization scale
  //   else printf( "I dont know what to do.\n");
  
  double mu = 0;

  if		( Func == kScale1 )		mu	= scale1 ;
  else if	( Func == kScale2 )		mu	= scale2 ;
  else if	( Func == kQuadraticSum )	mu	= FuncMixedOver1(scale1,scale2) ;
  else if	( Func == kQuadraticMean )	mu	= FuncMixedOver2(scale1,scale2) ;
  else if	( Func == kQuadraticSumOver4 )	mu	= FuncMixedOver4(scale1,scale2) ;
  else if	( Func == kScaleMax )		mu	= FuncMax(scale1,scale2);
  else if	( Func == kScaleMin )		mu	= FuncMin(scale1,scale2);
  else printf( "Error. could not identify functional form for scales calculation.\n");
  
  return scalefac * mu;

}


//______________________________________________________________________________
double FastNLOReader::FuncMixedOver1(double scale1 , double scale2 ){
  return ( sqrt( (pow(scale1,2) + pow(scale2,2))  / 1. ) ) ;
}

//______________________________________________________________________________
double FastNLOReader::FuncMixedOver2(double scale1 , double scale2 ){
  return ( sqrt( (pow(scale1,2) + pow(scale2,2))  / 2. ) ) ;
}

//______________________________________________________________________________
double FastNLOReader::FuncMixedOver4(double scale1 , double scale2 ){
  return ( sqrt( (pow(scale1,2) + pow(scale2,2))  / 4. ) ) ;
}

//______________________________________________________________________________
double FastNLOReader::FuncLinearMean(double scale1 , double scale2 ){
  return ( scale1 + scale2 ) / 2. ;
}

//______________________________________________________________________________
double FastNLOReader::FuncLinearSum(double scale1 , double scale2 ){
  return scale1 + scale2;
}

//______________________________________________________________________________
double FastNLOReader::FuncMax(double scale1 , double scale2 ){
  if ( scale1 > scale2 ) return scale1;
  else return scale2;
}

//______________________________________________________________________________
double FastNLOReader::FuncMin(double scale1 , double scale2 ){
  if ( scale1 < scale2 ) return scale1;
  else return scale2;
}




//______________________________________________________________________________



double FastNLOReader::SetScaleVariation(int scalevar , bool ReFillCache ){ 
  
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
    printf("    Please use SetScaleFacMuR(double) and SetScaleFacMuF(double).\n");
  }
  

  if (  scalevar >= BlockB_NLO->Nscalevar[0]  ){
    printf("Warning in FastNLOReader::SetScaleVariation. This table has only %d scalevariations stored. You wanted to acces number %d. Using '0' instead.\n", BlockB_NLO->Nscalevar[0] ,scalevar );
    fScalevar	= 0;
    return BlockB_NLO->ScaleFac[0][0] ;
  }
  
  fScalevar	= scalevar;
  printf(" * FastNLOReader::SetScaleVariation. Scalefactor of %4.2f for the nominal scale is chosen.\n",BlockB_NLO->ScaleFac[0][fScalevar]);

  if ( ReFillCache ){
    FillAlphasCache();
    FillPDFCache();
  }

  return BlockB_NLO->ScaleFac[0][fScalevar];
  
}




//______________________________________________________________________________



void FastNLOReader::SetFunctionalForm( EScaleFunctionalForm func , FastNLOReader::EMuX kMuX  ){
  //
  //  For MuVar tables this method sets the functional form of
  //  the renormalization or the factorization scale.
  //     func:  Choose a pre-defined function
  //     kMuX:  is it for mu_r or for mu_f
  //

  if ( BlockB_NLO->NScaleDep != 3 ) {
    printf("FastNLOReader::SetFunctionalForm. Warning. This is not a MuVar table.\n");
    printf("      SetFunctionalForm has no impact.\n");
    printf("      Please use another file, if you want to change your scale-definition.\n");
  }

  if ( kMuX == kMuR ) fMuRFunc = func;
  else fMuFFunc = func;

  if	( func == kScale2 || func == kQuadraticSum||  func == kQuadraticMean ||  func == kQuadraticSumOver4 ||  func == kScaleMax|| func == kScaleMin ) {
    if ( BlockB_LO->NscalenodeScalePt <= 3){
      printf("FastNLOReader::SetFunctionalForm. Error. There is no second scale variable available in this table.\n");
      printf("      Please use FastNLOReader::kScale1 only.\n");
      if ( kMuX == kMuR ) fMuRFunc = kScale1;
      else fMuFFunc = kScale1;
    }
    else if ( BlockB_LO->NscalenodeScalePt < 8 ){
      printf("FastNLOReader::SetFunctionalForm. Warning. Scale2 has only very little nodes (n=%d).\n",BlockB_LO->NscalenodeScalePt);
    }
  }
}


//______________________________________________________________________________


void FastNLOReader::SetMuRFunctionalForm( EScaleFunctionalForm func , bool ReFillCache  ){
  SetFunctionalForm(func,kMuR);

  if ( ReFillCache ){
    FillAlphasCache();
    FillPDFCache();
  }

}


//______________________________________________________________________________


void FastNLOReader::SetMuFFunctionalForm( EScaleFunctionalForm func , bool ReFillCache  ){
  SetFunctionalForm(func,kMuF);

  if ( ReFillCache ){
    FillPDFCache();
  }

}

//______________________________________________________________________________



void FastNLOReader::SetScaleFactorMuR(double fac , bool ReFillCache ){
  // 
  // Set scale factor for scale variations in MuVar tables
  // You have to ReFill your cache!
  // This is done automatically, but if you want to do it by yourself
  // set ReFillCache=false
  //

  if ( BlockB_NLO->NScaleDep != 3 ) {
    printf("FastNLOReader::SetScaleFactorMuR. Warning. This is not a MuVar table.\n");
    printf("      SetScaleFactorMuR has no impact.\n");
    printf("      Please use SetScaleVariation(int) instead.\n");
  }
  fScaleFacMuR = fac;

  if ( ReFillCache ){
    FillAlphasCache();
    FillPDFCache();
  }
}


//______________________________________________________________________________


void FastNLOReader::SetScaleFactorMuF(double fac , bool ReFillCache ){
  // 
  // Set scale factor for scale variations in MuVar tables
  // You have to ReFill your cache.
  // This is done automatically, but if you want to do it by yourself
  // set ReFillCache=false
  //
  
  if ( BlockB_NLO->NScaleDep != 3 ) {
    printf("FastNLOReader::SetScaleFactorMuF. Warning. This is not a MuVar table.\n");
    printf("      SetScaleFactorMuF has no impact.\n");
    printf("      Please use SetScaleVariation(int) instead.\n");
  }
  fScaleFacMuF = fac;

  if ( ReFillCache ){
    FillPDFCache();
  }
}


//______________________________________________________________________________



void FastNLOReader::ReadTable(void)
{
  //
  // Read in the FastNLO Table
  //
  
  // ---- check if file exists ----- //
  FILE* fp = fopen(ffilename.c_str(), "r");
  if (fp) {
    fclose(fp);
  } else {
    printf("Error. FastNLO table file does not exists. Was looking for: %s. Exiting.\n",ffilename.c_str());
    exit(1);
  } 
  
  // open stream
  ifstream* instream = new ifstream(ffilename.c_str(),ios::in);

  ReadBlockA1(instream);
  ReadBlockA2(instream);
  int nblocks	= Ncontrib + Ndata;
  //printf(" * FastNLOReader::ReadTable(). Reading %d B-Blocks.\n",nblocks);
  for(int i=0;i<nblocks;i++){
    FastNLOBlockB* blockb	= new FastNLOBlockB( "ReadingBlockB", NObsBin , instream );
    if ( blockb->IContrFlag1==1 && blockb->IContrFlag2==1 && blockb->IRef==0 ){
      if ( blockb->NScaleDep == 0 )   blockb->SetName("BlockB. LO. v2.0.");
      if ( blockb->NScaleDep == 3 )   blockb->SetName("BlockB. LO. v2.0 MuVar.");
      BlockB_LO		= blockb;
    }
    else if ( blockb->IContrFlag1==1 && blockb->IContrFlag2==2 && blockb->IRef==0 ){
      if ( blockb->NScaleDep == 0 )   blockb->SetName("BlockB. NLO. v2.0.");
      if ( blockb->NScaleDep == 3 )   blockb->SetName("BlockB. NLO. v2.0 MuVar.");
      BlockB_NLO	= blockb;
    }
    else if ( blockb->IContrFlag1==1 && blockb->IContrFlag2==1 && blockb->IRef==1 ){
      if ( blockb->NScaleDep == 0 )   blockb->SetName("BlockB. LO Reference. v2.0.");
      if ( blockb->NScaleDep == 3 )   blockb->SetName("BlockB. LO Reference. v2.0 MuVar.");
      BlockB_LO_Ref		= blockb;
    }
    else if ( blockb->IContrFlag1==1 && blockb->IContrFlag2==2 && blockb->IRef==1 ){
      if ( blockb->NScaleDep == 0 )   blockb->SetName("BlockB. NLO Reference. v2.0.");
      if ( blockb->NScaleDep == 3 )   blockb->SetName("BlockB. NLO Reference. v2.0 MuVar.");
      BlockB_NLO_Ref	= blockb;
    }
    else {
      printf("FastNLOReader::ReadTable(). Error in initializing the 'B'-Blocks. Your FastNLO Table file might be incompatible with this reader. Exiting.\n");
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
  //
  //  Read in information which is called Block A1 from file 
  //

  table->peek();
  if (table->eof()){
    printf("FastNLOReader::Read: Cannot read from file.\n");
    return;
  }
   
  int key = 0;
  *table >> key;
  if(key != tablemagicno){
    printf("FastNLOReader::Read: At beginning of block found %d instead of %d.\n",key,tablemagicno);
    return;
  };
  *table >> Itabversion;
  if ( Itabversion < 20000 ){
    printf("fnloBlockA1::Read. ERROR. This reader is only compatible with FastNLO v2.0 tables and higher.\n");  
    printf("       This FastNLO-table (file) is of version %6.4f.\n",Itabversion/10000.);
    printf("       Please download a compatible reader from the website or use the APPL_grid interface.\n");
    printf("       Exiting.\n");
    exit(1);
  }
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
    printf("FastNLOReader::Read: At end of block found %d instead of %d.\n",key,tablemagicno);
    return;
  };
  // Put magic number back
  for(int i=0;i<(int)(log10((double)key)+1);i++){
    table->unget();
  }
}


//______________________________________________________________________________


void FastNLOReader::ReadBlockA2(istream *table){
  //
  //  Read in information which is called Block A2 from file 
  //
  
  table->peek();
   if (table->eof()){
      printf("FastNLOReader::Read: Cannot read from file.\n");
      return;
   }

   int key = 0;
   *table >> key;
   if(key != tablemagicno){
      printf("FastNLOReader::Read: At beginning of block found %d instead of %d.\n",key,tablemagicno);
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
	  //irap++;
	  //cout << "irap: " << irap << ", RapIndex: " << RapIndex[irap] << endl;
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
      printf("FastNLOReader::Read: At end of block found %d instead of %d.\n",key,tablemagicno);
      return;
   };
   // Put magic number back
   for(int i=0;i<(int)(log10((double)key)+1);i++){
      table->unget();
   }
}



//______________________________________________________________________________


void FastNLOReader::Print(){
  //
  // Print FastNLO(Reader) internal variables
  //

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


void FastNLOReader::PrintCrossSections( ){
  //
  // Print Cross sections in NLO, k-factors and Reference table cross sections
  //
  
  if ( XSection.empty() && XSectionMuVar.empty() ){
     CalcCrossSection();
  }
  if ( XSectionRef.empty() && XSectionRefQ2.empty() ){
    CalcReferenceCrossSection();
  }

  vector < double > xs;
  vector < double > xsref;

  if ( BlockB_NLO->NScaleDep == 3 ) xs = XSectionMuVar;
  else xs = XSection;
  
  if ( BlockB_NLO->NScaleDep == 3 ){
    if ( fMuFFunc == kScale1 && fMuRFunc == kScale1 )			xsref = XSectionRefQ2;
    else if ( fMuFFunc == kScale1 && fMuRFunc == kQuadraticMean )	xsref = XSectionRefMufQ2MuRMixed;
    else if ( fMuFFunc == kQuadraticMean && fMuRFunc == kQuadraticMean )xsref = XSectionRefMixed;
    else xsref = XSectionRefMixed;
  }
  else xsref = XSectionRef;


  printf(" *  \n");
  printf(" *  FastNLO Cross sections for\n");
  for ( int i = 0 ; i < NScDescript ; i++ ){
    printf(" *     %s\n",ScDescript[i].c_str());
  }
  printf(" *  at sqrt(s) = %8.2f GeV\n", Ecms);
  printf(" *  \n");
  printf(" *  This is a %s-differential table in %s", ( (NDim==1)?"single":"double"),DimLabel[0].c_str());
  if ( NDim==2 ) printf(" and %s",DimLabel[1].c_str());
  printf(".\n");
  printf(" *\n");

  string unit[16] = { "b","","","mb","","","mu b","","","nb","","","pb","","","fb"};

  if ( NDim == 2 ){
    double lobindim2 = -321312;
    printf(" *  - Bin - |   ---  %s  ---         -- XS-FNLO [%s] --  -- k-factor -- |  -- XS-ref (NLOJET++) --    Diff [%%]\n",DimLabel[0].c_str(),unit[Ipublunits].c_str());
    printf(" *  ------------------------------------------------------------------------------------------------------------\n");
    for ( unsigned int i=0;i<xs.size();i++){
      if ( LoBin[i][1] != lobindim2 ){
	printf(" *                    ---->  from %9.3f to %9.3f in %s  <----\n",LoBin[i][1],UpBin[i][1],DimLabel[1].c_str());
	lobindim2 = LoBin[i][1];
      }
      printf(" *   %4.0f   | %9.3f - %9.3f     %12.4f          %5.3f      |     %12.4f            %5.4f\n",i*1.,LoBin[i][0],UpBin[i][0],xs[i],kFactor[i],xsref[i],(xs[i]-xsref[i])/xsref[i]*100.);
    }
  }

  else {
    printf("FastNLOReader::PrintCrossSections( ). Info. Single differential printing of cross sections not yet nicely implemented.\n");
    printf("   ---  %s  ---        - Bin -    -- XS-FNLO [pb] --    -- XS-ref (NLOJET++) --    Diff [%%]\n",DimLabel[0].c_str());
    for ( unsigned int i=0;i<xs.size();i++){
      printf("  %9.3f - %9.3f   %3.0f         %12.4f           %12.4f          %5.4f\n",LoBin[i][0],UpBin[i][0],i*1.,xs[i],xsref[i],(xs[i]-xsref[i])/xsref[i]*100.);
    }
  }
  printf(" *  ------------------------------------------------------------------------------------------------------------\n");
}


//______________________________________________________________________________


vector < double > FastNLOReader::GetXSection( ){
  // Get fast calculated NLO cross section
  
  if ( XSection.empty() && XSectionMuVar.empty() ){
    CalcCrossSection();
  }
  
  if ( BlockB_NLO->NScaleDep == 3 )
    return XSectionMuVar;
  else return XSection;
}


//______________________________________________________________________________


vector < double > FastNLOReader::GetReferenceXSection( ){
  // Get reference cross section from direct nlojet++ calculation
  
  if ( XSectionRef.empty() && XSectionRefQ2.empty() ){
    CalcReferenceCrossSection();
  }
  
  if ( BlockB_NLO->NScaleDep == 3 ){
    if ( fMuFFunc == kScale1 && fMuRFunc == kScale1 )			return XSectionRefQ2;
    else if ( fMuFFunc == kScale1 && fMuRFunc == kQuadraticMean )	return XSectionRefMufQ2MuRMixed;
    else if ( fMuFFunc == kQuadraticMean && fMuRFunc == kQuadraticMean )return XSectionRefMixed;
    else return XSectionRefMixed;
  }
  else return XSectionRef; // XSectionRef from BlockB-Ref
    
}


//______________________________________________________________________________


void FastNLOReader::CalcReferenceCrossSection( ){
  //
  //  Initialize the internal arrays for the reference cross
  //  sections with the information from the FastNLO file
  //
  
  XSectionRef.clear();
  XSectionRef.resize(NObsBin);
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
  
  if ( BlockB_LO->NScaleDep != 3 && ( BlockB_NLO_Ref==NULL ) )
    printf("FastNLOReader::CalcReferenceXSection( ). Warning. No reference cross sections available.\n");
       
}


//______________________________________________________________________________


void FastNLOReader::CalcCrossSection( ){
  //
  //  Initialize the internal arrays with the NLO cross sections
  //  with the information from the FastNLO file, the pdf and
  //  the defined alpha_s
  //
  
  XSection_LO.clear();
  XSection.clear();
  XSectionMuVar_LO.clear();
  XSectionMuVar.clear();
  kFactor.clear();

  XSection_LO.resize(NObsBin);
  XSection.resize(NObsBin);
  if ( BlockB_LO->NScaleDep == 3 ){
    XSectionMuVar_LO.resize(NObsBin);
    XSectionMuVar.resize(NObsBin);
  }
  kFactor.resize(NObsBin);
  
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

    if ( BlockB_LO->NScaleDep == 3 ){
      
      for(int jQ=0;jQ<BlockB_LO->NscalenodeScaleQ;jQ++){
	for(int jPt=0;jPt<BlockB_LO->NscalenodeScalePt;jPt++){
		    
	  double Q2   = BlockB_LO->ScaleNodeQ[i][jQ]*BlockB_LO->ScaleNodeQ[i][jQ];
	  
	  double mur	= CalcMu( kMuR , BlockB_LO->ScaleNodeQ[i][jQ] ,  BlockB_LO->ScaleNodePt[i][jPt] , fScaleFacMuR );
	  double muf	= CalcMu( kMuF , BlockB_LO->ScaleNodeQ[i][jQ] ,  BlockB_LO->ScaleNodePt[i][jPt] , fScaleFacMuF );

	  double mur2 = mur*mur;
	  double muf2 = muf*muf;

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
		    }
	}
      }
    }
    
  }
  
  // ---- k-factor calculation ---- //
  for(int i=0;i<NObsBin;i++){
    if ( BlockB_LO->NScaleDep == 3 )		kFactor[i]	= XSectionMuVar[i] / XSectionMuVar_LO[i];
    else					kFactor[i]	= XSection[i] / XSection_LO[i];
  }

}

//______________________________________________________________________________


void FastNLOReader::SetAlphasMz( double AlphasMz , bool ReCalcCrossSection ){
  //
  //  Set the alpha_s value at M_Z
  //
  
  if ( AlphasMz != fAlphasMz ){
    fAlphasMz	= AlphasMz;		// new alpha_s value
    FillAlphasCache();
    if ( ReCalcCrossSection ) CalcCrossSection(); 
  }
  else {
    // nothing to do!
  }
  
}



//______________________________________________________________________________


void FastNLOReader::FillAlphasCache(){
  //
  //  Fill the internal alpha_s cache.
  //  This is usally called automatically. Only if you
  //  make use of ReFillCache==false options, you have
  //  to take care of this filling by yourself.
  //

  FillAlphasCacheInBlockB( BlockB_LO  );
  FillAlphasCacheInBlockB( BlockB_NLO );

}



//______________________________________________________________________________


double FastNLOReader::GetAlphas( double Q ){
  // 
  //  Internal method for caluclating the alpha_s(mu)
  //
  
  //switch (AlphasEvolution )
  if ( fAlphasEvolution == kGRV )			return GetAlphasGRV	( Q , fAlphasMz );
  else if ( fAlphasEvolution == kNLOJET )		return GetAlphasNLOJET	( Q , fAlphasMz );
  else if ( fAlphasEvolution == kCTEQpdf )		return GetAlphasCTEQpdf	( Q , fAlphasMz );
  else if ( fAlphasEvolution == kFastNLO )		return GetAlphasFastNLO	( Q , fAlphasMz );
  else if ( fAlphasEvolution == kLHAPDFInternal )	return GetAlphasLHAPDF	( Q );
  else if ( fAlphasEvolution == kQCDNUMInternal )	return GetAlphasQCDNUM	( Q );
  else return 0;
}



//______________________________________________________________________________


void FastNLOReader::FillAlphasCacheInBlockB( FastNLOBlockB* B ){
  // 
  //  Internal method for filling alpha_s cache
  //
  

  for(int i=0;i<NObsBin;i++){
    for(int j=0;j<B->GetTotalScalenodes();j++){
      int scalenode1 = j;
      int scalenode2 = j;
      if (B->NScaleDim>1){          
	scalenode1 = j / B->Nscalenode[1];
	     scalenode2 = j % B->Nscalenode[1];
      }
      
      double mur	= B->ScaleNode[i][0][fScalevar][scalenode1];
      double as		= GetAlphas(mur);
      
      double alphastwopi = pow( as/TWOPI , B->Npow );
      B->AlphasTwoPi_v20[i][j] = alphastwopi;
    }
  
    if ( B->NScaleDep == 3 ){
	for(int jQ=0;jQ<B->NscalenodeScaleQ;jQ++){
	  for(int jPt=0;jPt<B->NscalenodeScalePt;jPt++){
	    // 	    double Q2   = B->ScaleNodeQ[i][jQ]*B->ScaleNodeQ[i][jQ];
	    // 	    double Pt   = B->ScaleNodePt[i][jPt];
	    double mur		= CalcMu( kMuR , BlockB_LO->ScaleNodeQ[i][jQ] ,  BlockB_LO->ScaleNodePt[i][jPt] , fScaleFacMuR );
	    double as		= GetAlphas(mur);
	    double alphastwopi	= pow( as/TWOPI, B->Npow );
	    B->AlphasTwoPi[i][jQ][jPt] = alphastwopi;
	  }
	}
    }
  }
}


//______________________________________________________________________________


double FastNLOReader::GetAlphasLHAPDF(double Q){
  //
  // Implementation of Alpha_s evolution as function of Mu_r only.
  //
  // the alpha_s evolution is done within LHAPDF.
  // 
  // WARNING: You cannot change alpha_s(Mz), but is is
  // defined with the pdf.
  //
  
  return LHAPDF::alphasPDF(Q); 
}


//______________________________________________________________________________


double FastNLOReader::GetAlphasQCDNUM(double Q){
   //
   // Sorry! This function is still under developement.
   printf("FastNLOReader::GetAlphasQCDNUM. ERROR. Not yet implemented.\n");
   return 0.;
}


//______________________________________________________________________________


double FastNLOReader::GetAlphasNLOJET(double Q, double alphasMZ){
  //
  // Implementation of Alpha_s evolution as function of Mu_r.
  //
  // this is the evolution, which is used by nlojet++ and cteq6m.
  // Be aware of the Mz-value from 2001 of 91.70.
  //
  // Values used in nlojet++ and cteq6m:
  //   alphasmz   = 0.1179  
  //   double b0  = 1.2202
  //   double b1  = 0.4897
  //   double Mz  = 91.70
  //
  // as evolution by lhpdf.c 
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
  double L = log(Q/Mz);
  L = (b0 + alphasMZ*b1)*L;

  return alphasMZ/(1.0 + alphasMZ*L);

}


//______________________________________________________________________________


double FastNLOReader::GetAlphasGRV(double MU, double ALPSMZ){
  //
  // Implementation of Alpha_s evolution as function of Mu_r.
  //
  // *******************************************************************
  //    the version from v1.4 as used in fnh2003-desy07073 HERA-I HighQ2 jets
  //    from: fn-alphas-demo.f 
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
  // *******************************************************************
  // * This c++ translation:
  // * D. Britzger  03/03/11
  // *

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

}


//______________________________________________________________________________


double FastNLOReader::GetAlphasCTEQpdf(double Q, double alphasMZ){
  //
  // Implementation of Alpha_s evolution as function of Mu_r.
  //
  // alpha_s evolution as it is within by cteq-pdf-1.0.4 and used in nlojet 4.1.3 without crack

  double as_twopi = alphasMZ/TWOPI;
  double Mz	= 91.187;
  //double Mz	= 91.70;
  //   double Mz	= MZ;

  //int ord=2;
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
  // Implementation of Alpha_s evolution as function of Mu_r.
  //
  // This is the 'original' FNLOv2.0 implementation.
  //

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


void FastNLOReader::FillPDFCache( bool ReCalcCrossSection ){
  //
  //  Fill the internal pdf cache.
  //  This function has to be called by the user, since the 
  //  pdf parameters and evolutions are calculated externally.
  //
   
  if ( fPDFInterface == kLHAPDF ){
     if ( fLHAPDFfilename == "" || fLHAPDFpath =="" ){
	printf("FastNLOReader::FillPDFCache(). ERROR. You must specify a LHAPDF filename first or you have to specify kH1FITTER..\n"); exit(1);
     }
     InitLHAPDF();
  }
  else if ( fPDFInterface == kH1FITTER ){
     evolution_();
  }
  
  FillBlockBPDFLCs(BlockB_LO);
  FillBlockBPDFLCs(BlockB_NLO);
 
  if ( ReCalcCrossSection ) CalcCrossSection();
 
}



//______________________________________________________________________________


void FastNLOReader::InitLHAPDF(){
  //
  //  Initalize some necessary LHAPDF parameters
  //
    
  if ( fLHAPDFfilename == "" || fLHAPDFpath ==""){
    printf("FastNLOReader::FillPDFCacheLHAPDF(). ERROR. You must specify a LHAPDF filename first.\n"); exit(1);
  }

  string LHAPDFfile = fLHAPDFpath+"/"+fLHAPDFfilename;
  
  // ---- check if file exists ----- //
  FILE* fp = fopen(LHAPDFfile.c_str(), "r");
  if (fp) {
    fclose(fp);
  } else {
    printf("Error. LHAPDF file does not exists. Was looking for: %s. Exiting.\n",LHAPDFfile.c_str());
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


void FastNLOReader::FillBlockBPDFLCs( FastNLOBlockB* B ){
  
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
	   xfx = GetXFX(xp,muf);
	   //xfx = LHAPDF::xfx(xp,muf); // LHAPDF::xfx_p_(x,muf,0,0)
	   vector < double > buffer = CalcPDFLinearComb(xfx,xfx,B->IPDFdef1, B->IPDFdef2, B->NSubproc ); //calculate linear combinations
	   for(int l=0;l<B->NSubproc;l++){ 
	     B->PdfLc[i][j][k][l] = buffer[l];
	   }
	 }
       }
     }
		
     if ( B->NScaleDep == 3 ){
		xfx.clear();// do i need this?
		xfx.resize(13);// do i need this?
		for(int i=0;i<NObsBin;i++){
		  int nxmax = B->GetNxmax(i);
		  for(int jQ=0;jQ<B->NscalenodeScaleQ;jQ++){
		    for(int jPt=0;jPt<B->NscalenodeScalePt;jPt++){
		      double muf	= CalcMu( kMuF , BlockB_LO->ScaleNodeQ[i][jQ] ,  BlockB_LO->ScaleNodePt[i][jPt] , fScaleFacMuF );
		      
		      for(int x=0;x<nxmax;x++){ 
			double xp	= B->XNode1[i][x];
			xfx = GetXFX(xp,muf);
			//xfx = LHAPDF::xfx(xp, muf); // LHAPDF::xfx_p_(x,muf,0,0)
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

   else {
      printf("FastNLOReader::FillBlockBPDFLCs. Error. I don't know what to do.\n");
   }
   


}



//______________________________________________________________________________



vector<double> FastNLOReader::GetXFX(double xp, double muf){
   //
   //  Internal method.
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   // 
  
   if ( fPDFInterface == kLHAPDF ){
      return LHAPDF::xfx(xp,muf);
   }
   else if ( fPDFInterface == kH1FITTER ){
      int iqnset = 1;
      int iqnchk = 0;
      double muf2	= muf*muf;
      vector < double > a;
      a.resize(13);
      fpdfxq_(&iqnset, &xp, &muf2, &a[0], &iqnchk); 
      return a;
   }
   else {
      vector < double > a;
      return a;
   }
}



//______________________________________________________________________________


vector<double> FastNLOReader::CalcPDFLinearComb(vector<double> pdfx1, vector<double> pdfx2, int IPDFdef1, int IPDFdef2, int NSubproc){
   //
   //  Internal method.
   //  CalcPDFLinearComb is used to calculate the 
   //  linear combinations of different parton flavors
   //  according to the used subprocesses of your
   //  calculation/table.
   //
   
   vector < double > ret;

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
	 printf("FastNLOReader::CalcPDFLinearComb. Error no support for pp and ppbar tables.\n");
	 break;
      default: printf("fnloBlockB::CalcPDFLinearComb ::Ipdfdef1=4, Ipdfdef2= %d not supported. Exit.\n",IPDFdef2); exit(1); return ret;
      }
      break;
   default: printf("fnloBlockB::CalcPDFLinearComb :Ipdfdef1= %d not supported. Exit.\n",IPDFdef1); exit(1); return ret;
   }

  
}


//______________________________________________________________________________


vector<double> FastNLOReader::CalcPDFLinearCombDIS( vector<double> pdfx1 , int NSubproc){
   //
   //  Internal method.
   //  CalcPDFLinearComb is used to calculate the 
   //  linear combinations of different parton flavors
   //  according to the used subprocesses of your
   //  calculation/table.
   //
  
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


