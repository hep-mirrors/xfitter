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
  Int_t Prepare();
  void PrintCanvas(TCanvas* can);
  Int_t PlotPdfSub(TVirtualPad* pad, H1FitterOutput* FitterOut,H1FitterOutput* FitterRef,
		   Int_t Q2Bin, const Char_t* Title, H1FitterOutput::pdf pdf1, H1FitterOutput::pdf pdf2,
		   TVirtualPad* legend, TObjArray* TrashBin);

  void ScaleGraph2ToGraph1(TGraph* graph1, TGraph* graph2, TLine*& line, TGaxis*& axis, Double_t MeanRatio);
 public:
  H1FitterPainter(Int_t  Bands = 0);
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

  // Number of PDF eighenvectors for PDF bands.
  Int_t nBands;
};

#endif
