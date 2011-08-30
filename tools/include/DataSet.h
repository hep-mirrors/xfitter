#ifndef  DataSet_h
#define  DataSet_h

#include <TString.h>
#include <TObject.h>
#include <iostream>
#include <TGraphErrors.h>
#include <TObjArray.h>
#include <vector>
#include <TObjString.h>
#include <TOrdCollection.h>

using std::cout;
using std::cerr;
using std::endl;


class  DataSet {
 private:    
  TString* fName;
  Int_t fSetId;
  TString* fV1;
  TString* fV2;
  TString* fV3;
  TOrdCollection* fDUnc;
  TOrdCollection* fDTot;
  TOrdCollection* fTheo;
  TOrdCollection* fLabels;
  std::vector <Double_t> fDLabels;

 public: 
  DataSet();
  DataSet(Int_t SetId, const Char_t* name, const Char_t* v1, const Char_t* v2, const Char_t* v3);
  virtual ~DataSet();
  Int_t inline GetSetId() {return fSetId;};
  TGraphErrors* GetDataUncr(Int_t i);
  TGraphErrors* GetDataTotal(Int_t i);
  TGraphErrors* GetTheory(Int_t i);
  TObjString*   GetLabel(Int_t i);
  void AddPoint(Double_t v1, Double_t v2, Double_t v3, Double_t data, Double_t uncorrerr, Double_t toterr, Double_t theory, Double_t pull);
  Double_t GetMinimum(Int_t i);
  Double_t GetChi2(Int_t i);
  Int_t GetNpts(Int_t i);
  inline const Char_t* GetName() {return fName->Data();}
  Int_t FillLabels(TH1F* h);
  inline Int_t GetNGraphs() {return fDUnc->GetEntries();}
   
 private:
  Int_t FindGraphIndex(Double_t value, const Char_t* label);
  void AddPoint(Int_t GraphIdx, Double_t x, Double_t data, Double_t uncorrerr, Double_t toterr, Double_t theory);

      };
#endif
