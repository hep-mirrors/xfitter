#ifndef  H1FitterOutput_h
#define  H1FitterOutput_h

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

using std::cout;
using std::cerr;
using std::endl;

const Int_t NColumn = 14;

class  H1FitterOutput {
 public:
  enum pdf{kNULL=-1, kGluon=0, kU=1, kD=2, kUv=3, kDv=4, kUb=5, kDb=6, kSea=7, kS=8, kC=9, kB=10};

 private:    
  TString* fDirectory;
  static const Int_t fNpdfs = 11;   
  static const Int_t fNBands = 20;   
  Int_t fNpoints;                   // Number of x-points 
  Int_t fNQ2Files;                  // Number of Q2 files
  TObjArray* fPdfs[fNpdfs];
  Double_t   fQ2Value[100];         // up to 100 files
  TH1F*      fPull;
  
  Int_t fNDataSets;
  std::vector<DataSet*> fDataSets;

 public:
   H1FitterOutput(const Char_t* directory);
   virtual ~H1FitterOutput();
   Int_t Prepare(bool DrawBand);
   TGraph* GetPdf(H1FitterOutput::pdf ipdf, Int_t Q2bin);
   inline TString* GetDirectory() {return fDirectory;}
   Int_t GetNsets();
   DataSet* GetSet(Int_t i) {return fDataSets[i];}
   inline TH1F* GetPull() {return fPull;}
   
   // Return number of Q2 files
   const Int_t GetNQ2Files() { return fNQ2Files;}
   // Return Q2 value for a given file
   const Double_t GetQ2Value(Int_t iQ2bin);

 private:
   Int_t PrepareDataSets();
   Int_t PreparePdf(bool DrawBand);
   void SetPdfPoint(Int_t ipdf, Int_t iq2, Int_t ipoint, Double_t x, Double_t y);
   void SetPdfError(Int_t ipdf, Int_t iq2, Int_t ipoint, Double_t x, Double_t Up, Double_t Dn);
   Bool_t CheckDirectory();
};
#endif
