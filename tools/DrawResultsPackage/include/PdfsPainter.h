#ifndef PdfsPainter_h
#define PdfsPainter_h

#include <vector>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>

using namespace std;

  struct gstruct{
    TGraphAsymmErrors* graph;
    string label;
  };

TCanvas * PdfsPainter(double q2, string pdf, vector <gstruct> pdfgraphs);
TCanvas * PdfsRatioPainter(double q2, string pdf, vector <gstruct> pdfgraphs);

#endif
