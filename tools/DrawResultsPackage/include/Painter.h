#ifndef Painter_h
#define Painter_h


#include <stdlib.h>
#include <iostream>
#include <TString.h>
#include <Output.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLegend.h>
#include <TGaxis.h>
#include <DataSet.h>
#include <Rtypes.h>
#include <TObjArray.h>
#include <TPaveText.h>

using std::cout;
using std::cerr;
using std::endl;
using std::flush;

class Painter  {

 protected:
  Painter(Painter &foo);
  Painter& operator=(Painter &foo);
  Int_t DrawPDF(Int_t ival);
  Int_t DrawPull();
  
  Int_t DrawDataSet(DataSet* dataset, DataSet* datasetref, Bool_t RatioToData, EColor=kRed);
  Int_t DrawDataSetEMP(DataSet* datasetB, DataSet* datasetE, DataSet* datasetM, DataSet* datasetP, Bool_t RatioToData);
  
  Int_t DrawDataSetRatio(DataSet* dataset, DataSet* datasetref);
  
  Int_t DrawFitResults();
  Int_t DrawFitResultsEMP();
  
  void DrawMessages(Output* output);
  void DrawCorrelations(Output* output);
  Int_t Prepare();
  void PrintCanvas(TCanvas* can);
  //TCanvas* PrepareDataSetCanvas(DataSet* dataset, DataSet* datasetref, TObjArray* TrashBin, Double_t& MarkerSize, Bool_t Ratio);
  
  Int_t PlotPdfSub(TVirtualPad* pad, Output* FitterOut,Output* FitterRef,
		   Int_t Q2Bin, const Char_t* Title, Output::pdf pdf1, Output::pdf pdf2,
		   TVirtualPad* legend, TObjArray* TrashBin);

  Int_t PlotPdfSubEMP(TVirtualPad* pad, Int_t Q2Bin, const Char_t* Title, Output::pdf pdf1, Output::pdf pdf2,
		   TVirtualPad* legend, TObjArray* TrashBin);
  
  void ScaleGraph2ToGraph1(TGraph* graph1, TGraph* graph2, TLine*& line, TGaxis*& axis, Double_t MeanRatio);
  TText* AddLineToPave(TObjArray* paves, float& yposition, const char* text, const char* option);
  void FillPavesWithFitResults(TObjArray* paves, Output* output);
 public:
  Painter(bool  Bands = false, bool DrawBase = false, bool DrawExp = false, bool DrawModel = false, bool DrawParam = false);
  virtual ~Painter();
  Int_t Draw();


  inline void SetPath(const Char_t* path) {fPath->Form(path);}
  inline void SetPathRef(const Char_t* path) {fPathRef->Form(path);}
  inline void SetPathBase(const Char_t* path) {fPathBase->Form(path);}
  inline void SetPathExp(const Char_t* path) {fPathExp->Form(path);}
  inline void SetPathModel(const Char_t* path) {fPathModel->Form(path);}
  inline void SetPathParam(const Char_t* path) {fPathParam->Form(path);}

 protected: 
  TString* fPath;
  TString* fPathRef;
  Output* fOutput;
  Output* fOutputRef;
  TString* fPsFileName;
  EColor fColor;
  EColor fColorRef;
  EColor fHighlight;
  Int_t  fFillStyle;
  Int_t  fFillStyleRef;
  
  bool fDrawBase;
  bool fDrawExp;
  bool fDrawModel;
  bool fDrawParam;
  TString* fPathBase;
  TString* fPathExp;
  TString* fPathModel;
  TString* fPathParam;
  Output* fOutputBase;
  Output* fOutputExp;
  Output* fOutputModel;
  Output* fOutputParam;
  EColor fColorBase;
  EColor fColorExp;
  EColor fColorModel;
  EColor fColorParam;
  Int_t  fFillStyleBase;
  Int_t  fFillStyleExp;
  Int_t  fFillStyleModel;
  Int_t  fFillStyleParam;
  

  // Number of PDF eighenvectors for PDF bands.
  bool fBands;
};

#endif
