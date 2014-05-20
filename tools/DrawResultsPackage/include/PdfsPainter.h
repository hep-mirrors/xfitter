#ifndef PdfsPainter_h
#define PdfsPainter_h

#include <vector>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>

#include "PdfData.h"

using namespace std;

extern vector <string> pdflabels;
extern vector <string> pdffiles;

vector <TCanvas*> PdfsPainter(double q2, pdftype ipdf);
//TCanvas * PdfsRatioPainter(double q2, pdftype ipdf);

#endif
