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


class  H1FitterOutput {
 public:
  enum pdf{kNULL=-1, kGluon=0, kU=1, kD=2, kUv=3, kDv=4, kUb=5, kDb=6, kSea=7, kS=8, kC=9, kB=10};

 private:    
  TString* fDirectory;
  static const Int_t fNpdfs=11;
  static const Int_t fNpoints=100;
  TObjArray* fPdfs[fNpdfs];
  TH1F*      fPull;
  
  Int_t fNDataSets;
  std::vector<DataSet*> fDataSets;

 public:
   H1FitterOutput(const Char_t* directory);
   virtual ~H1FitterOutput();
   Int_t Prepare();
   TGraph* GetPdf(H1FitterOutput::pdf ipdf, Int_t Q2bin);
   inline TString* GetDirectory() {return fDirectory;}
   Int_t GetNsets();
   DataSet* GetSet(Int_t i) {return fDataSets[i];}
   inline TH1F* GetPull() {return fPull;}
 private:
   Int_t PrepareDataSets();
   Int_t PreparePdf();
   void SetPdfPoint(Int_t ipdf, Int_t iq2, Int_t ipoint, Double_t x, Double_t y);
   Bool_t CheckDirectory();
};
#endif
