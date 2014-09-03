#ifndef PdfData_h
#define PdfData_h

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <TGraphAsymmErrors.h>

using namespace std;

enum pdferr {AsymHess, SymHess, MC};
enum pdftype{uv=0, dv, g, Sea, ubar, dbar, s, Rs, c, b, dbarubar, uvdv, U, D, Ubar, Dbar};

extern vector <pdftype> pdfs;
extern vector <string> pdflabels;
extern vector <string> pdffiles;

//Class storing one PDF table
class Pdf
{
 public:
  Pdf() {};
  Pdf(string filename);
  // Copy constructor: (not used)
  /*  Pdf(const Pdf& prior) : Q2value (prior.Q2value),
    Xmin (prior.Xmin),
    Xmax (prior.Xmax),
    NxValues ( prior.NxValues),
    NPdfs    ( prior.NPdfs),
    tablemap ( prior.tablemap),
    tablemapup ( prior.tablemapup),
    tablemapdn ( prior.tablemapdn),
    xpoints  (prior.xpoints)      {}*/

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

struct pdfshift 
{
  double val;
  double err;
  int    id;
};


//Class storing all the PDF data of one directory
class PdfData
{
 public:
  PdfData() {};
  PdfData(string dirname, string label);
  // Copy constructor, to for profiled PDF: (not used)
  //PdfData(const PdfData &Prior, string dirname, string label);
  void profile(string dirname, string label);  //profile PDF uncertainty bands

  pdferr err;   //Type of PDF uncertainties
  bool model;   //Model PDF uncertainties
  bool par;   //Parametrisation PDF uncertainties

  map <float, Pdf> Central;               //Q2 values central PDF map
  map <float, vector <Pdf> > Errors;      //Q2 values PDF errors map
  map <float, vector <Pdf> > ModelErrors; //Q2 values PDF errors map
  map <float, vector <Pdf> > ParErrors;   //Q2 values PDF errors map
  map <float, Pdf> Up;                    //Q2 values up PDF map
  map <float, Pdf> Down;                  //Q2 values down PDF map

  vector<pdfshift> pdfshifts;  // For profiled PDFs, keep shift info
};

#endif
