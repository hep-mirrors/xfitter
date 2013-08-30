#ifndef  Output_h
#define  Output_h

#include <TObject.h>
#include <iostream>
#include <TString.h>
#include <TObjArray.h>
#include <TGraph.h>
#include <DataSet.h>
#include <TArrayF.h>
#include <TArrayI.h>
#include <vector>
#include <TH1F.h>
#include <TGraphAsymmErrors.h>

using std::cout;
using std::cerr;
using std::endl;

const Int_t NColumn = 14;

class  Output {
 public:
  enum pdf{kNULL=-1, kGluon=0, kU=1, kD=2, kUv=3, kDv=4, kUb=5, kDb=6, kSea=7, kUSea=8, kDSea=9, kS=10, kC=11, kB=12};

 protected:    
  TString* fName;
  TString* fDirectory;
  static const Int_t fNpdfs = 13;
  static const Int_t fNBands = 20;   
  static const Int_t fMaxParameters = 50;
  Int_t fNpoints;                   // Number of x-points 
  Int_t fNQ2Files;                  // Number of Q2 files
  TObjArray* fPdfs[fNpdfs];
  Double_t   fQ2Value[100];         // up to 100 files
  TH1F*      fPull;
  
  Int_t fNDataSets;
  std::vector<DataSet*> fDataSets;

  TObjArray* fMessages;
  Double_t   fFittedParameters[fMaxParameters][2];
  TString*   fFittedParametersNames[fMaxParameters];
  Int_t      fNFittedParameters;
  Int_t      fNNuisanceParameters;
  Double_t   fCorrPar[fMaxParameters][fMaxParameters];

  TString*   fNuisanceParNames[fMaxParameters];
  Double_t   fNuisancePar[fMaxParameters][2];
  TString*   fErrorCalculationMethod;
  TString*   fCorrelationCalculationMethod;
  TString*   fErrorTrustLevel;
  Bool_t     fConverged;
  Bool_t     fFinished;
  Double_t   fChi2UncTotal;
  Double_t   fChi2CorTotal;

  Bool_t fParametersCheck;
  Bool_t fNuisanceCheck;
  Bool_t fCovarianceCheck;
  Bool_t fMessagesCheck;

 public:
   Output(const Char_t* directory);
   virtual ~Output();
   virtual Int_t Prepare(bool DrawBand, TString option = TString("b"));
   TGraphAsymmErrors* GetPdf(Output::pdf ipdf, Int_t Q2bin);
   inline TString* GetDirectory() {return fDirectory;}
   inline TString* GetName() {return fName;}
   Int_t GetNsets();
   DataSet* GetSet(Int_t i) {return fDataSets[i];}
   inline TH1F* GetPull() {return fPull;}
   inline Bool_t GetParametersCheck() {return fParametersCheck;}
   inline Bool_t GetNuisanceCheck() {return fNuisanceCheck;}
   inline Bool_t GetCovarianceCheck() {return fCovarianceCheck;}
   inline Bool_t GetMessagesCheck() {return fMessagesCheck;}
   inline TString* GetErrorTrustLevel() {return fErrorTrustLevel;}
   inline TString* GetErrorCalculationMethod() {return fErrorCalculationMethod;}
   inline TString* GetCorrelationCalculationMethod() {return fCorrelationCalculationMethod;}
   inline TObjArray* GetMessages() {return fMessages;}
   inline Double_t GetChi2UncTotal() {return fChi2UncTotal;}
   inline Double_t GetChi2CorTotal() {return fChi2CorTotal;}
   
   // Return number of Q2 files
   const Int_t GetNQ2Files() { return fNQ2Files;}
   // Return Q2 value for a given file
   const Double_t GetQ2Value(Int_t iQ2bin);
   inline TString* GetFittedParametersName(Int_t idx) {return fFittedParametersNames[idx];}
   Double_t GetFittedParameter(Int_t idx, Bool_t error=kFALSE);
   inline TString* GetNuisanceParNames(Int_t idx) {return fNuisanceParNames[idx];}
   inline Int_t GetNFittedParameters() {return fNFittedParameters;}
   inline Int_t GetNNuisanceParameters() {return fNNuisanceParameters;}
   Double_t GetNuisancePar(Int_t idx, Bool_t error=kFALSE);
   inline Double_t GetCorPar(int i, int j) { return fCorrPar[i][j];}
   Int_t PreparePdf(bool DrawBand, TString option = TString("b"));
   void PrepareParameters();
   inline Bool_t GetConverged() {return fConverged;}
   inline Bool_t GetFinished() {return fFinished;}

 protected:
   Int_t PrepareDataSets();
   void PrepareName();
   void SetPdfPoint(Int_t ipdf, Int_t iq2, Int_t ipoint, Double_t x, Double_t y);
   void SetPdfError(Int_t ipdf, Int_t iq2, Int_t ipoint, Double_t x, Double_t Up, Double_t Dn);
   Bool_t CheckDirectory();
   Bool_t CheckFile();
   void PrepareMandyParameters();
};
#endif
