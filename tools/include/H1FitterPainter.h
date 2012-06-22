#ifndef H1FitterPainter_h
#define H1FitterPainter_h


#include <stdlib.h>
#include <iostream>
#include <TString.h>
#include <H1FitterOutput.h>
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

class H1FitterPainter  {

 private:
  H1FitterPainter(H1FitterPainter &foo);
  H1FitterPainter& operator=(H1FitterPainter &foo);
  Int_t DrawPDF(Int_t ival);
  Int_t DrawPull();
  Int_t DrawDataSet(DataSet* dataset, DataSet* datasetref, EColor=kRed);
  Int_t DrawDataSetRatio(DataSet* dataset, DataSet* datasetref);
  Int_t DrawFitResults();
  void DrawMessages(H1FitterOutput* output);
  void DrawCorrelations(H1FitterOutput* output);
  Int_t Prepare();
  void PrintCanvas(TCanvas* can);
  Int_t PlotPdfSub(TVirtualPad* pad, H1FitterOutput* FitterOut,H1FitterOutput* FitterRef,
		   Int_t Q2Bin, const Char_t* Title, H1FitterOutput::pdf pdf1, H1FitterOutput::pdf pdf2,
		   TVirtualPad* legend, TObjArray* TrashBin);

  void ScaleGraph2ToGraph1(TGraph* graph1, TGraph* graph2, TLine*& line, TGaxis*& axis, Double_t MeanRatio);
  TText* AddLineToPave(TObjArray* paves, float& yposition, const char* text, const char* option);
  void FillPavesWithFitResults(TObjArray* paves, H1FitterOutput* output);
 public:
  H1FitterPainter(bool  Bands = false);
  virtual ~H1FitterPainter();
  Int_t Draw();


  inline void SetPath(const Char_t* path) {fPath->Form(path);}
  inline void SetPathRef(const Char_t* path) {fPathRef->Form(path);}

 private: 
  TString* fPath;
  TString* fPathRef;
  H1FitterOutput* fH1FitterOutput;
  H1FitterOutput* fH1FitterOutputRef;
  TString* fPsFileName;
  EColor fColor;
  EColor fColorRef;
  EColor fHighlight;
  Int_t  fFillStyle;
  Int_t  fFillStyleRef;
  

  // Number of PDF eighenvectors for PDF bands.
  bool fBands;
};

#endif
