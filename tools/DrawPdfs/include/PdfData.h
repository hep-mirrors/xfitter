#ifndef PdfData_h
#define PdfData_h

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <TGraphAsymmErrors.h>

using namespace std;

enum pdferr {AsymHess, SymHess, MC};
enum pdftype{g=0, uv, dv, ubar, dbar, s, U, D, Sea, Ubar, Dbar, c, b, dbarubar, sdbar};

extern vector <pdftype> pdfs;
extern vector <string> pdflabels;
extern vector <string> pdffiles;

//Class storing one PDF table
class Pdf
{
 public:
  Pdf() {};
  Pdf(string filename);

  int GetNx(){return NxValues;}     // Number of X points
  int GetNPdfs(){return NPdfs;}     // Number of PDfs
  double GetXmin(){return Xmin;}    // Xmin value
  double GetXmax(){return Xmax;}    // Xmax value
  float GetQ2() {return Q2value;}  // Q2

  TGraphAsymmErrors * GetPdf(pdftype ipdf);

  vector <double> GetTable(pdftype ipdf) {return tablemap[ipdf];};
  void SetPoint(pdftype ipdf, int ix, double value) {tablemap[ipdf][ix] = value;};
  void SetErrUp(pdftype ipdf, int ix, double value) {tablemapup[ipdf].resize(ix+1); tablemapup[ipdf][ix] = value;};
  void SetErrDn(pdftype ipdf, int ix, double value) {tablemapdn[ipdf].resize(ix+1); tablemapdn[ipdf][ix] = value;};

 private:
  float Q2value;
  double Xmin;
  double Xmax;
  int    NxValues;
  int    NPdfs;
  map <pdftype, vector <double> > tablemap;
  map <pdftype, vector <double> > tablemapup;
  map <pdftype, vector <double> > tablemapdn;
  vector <double> xpoints;
};

//Class storing all the PDF data of one directory
class PdfData
{
 public:
  PdfData() {};
  PdfData(string dirname, string label);

  pdferr err;   //Type of PDF uncertainties

  map <float, Pdf> Central;           //Q2 values central PDF map
  map <float, vector <Pdf> > Errors;  //Q2 values PDF errors map
  map <float, Pdf> Up;                //Q2 values up PDF map
  map <float, Pdf> Down;              //Q2 values down PDF map
};

#endif
