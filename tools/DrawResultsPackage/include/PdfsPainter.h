#ifndef PdfsPainter_h
#define PdfsPainter_h

#include <vector>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>

using namespace std;

extern vector <string> pdflabels;
extern vector <string> pdffiles;

  struct gstruct{
    TGraphAsymmErrors* graph;
    string label;
  };

TCanvas * PdfsPainter(double q2, int ipdf, vector <gstruct> pdfgraphs);
TCanvas * PdfsRatioPainter(double q2, int ipdf, vector <gstruct> pdfgraphs);

#endif
