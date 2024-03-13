#ifndef PdfData_h
#define PdfData_h

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <TGraphAsymmErrors.h>

using namespace std;

enum pdferr {None, AsymHess, SymHess, MC};
enum pdftype{uv=0, dv, g, Sea, ubar, dbar, s, sbar, Rs, soversbar, c, b, dbarminubar, uvmindv, U, D, Ubar, Dbar, goversea, doveru, dbaroverubar, dvoveruv, rs, photon, SeaOverGlue, photonOverGlue, uvplusdv, uvplusdvplusSea, uvplusdvNEW};

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

  TGraphAsymmErrors * GetPdf(pdftype ipdf, int iunctype = 3);
  TGraphAsymmErrors * GetPdfCen(pdftype ipdf) {return GetPdf(ipdf, 0);}
  TGraphAsymmErrors * GetPdfExp(pdftype ipdf) {return GetPdf(ipdf, 1);}
  TGraphAsymmErrors * GetPdfModel(pdftype ipdf) {return GetPdf(ipdf, 2);}
  TGraphAsymmErrors * GetPdfParam(pdftype ipdf) {return GetPdf(ipdf, 3);}

  vector <double> GetTable(pdftype ipdf) {return tablemap[ipdf];};
  void SetPoint(pdftype ipdf, int ix, double value) {tablemap[ipdf][ix] = value;};
  void SetErrUp(pdftype ipdf, int ix, double value) {tablemapup[ipdf].resize(ix+1); tablemapup[ipdf][ix] = value;};
  void SetErrDn(pdftype ipdf, int ix, double value) {tablemapdn[ipdf].resize(ix+1); tablemapdn[ipdf][ix] = value;};
  void SetExpUp(pdftype ipdf, int ix, double value) {tablemapexpup[ipdf].resize(ix+1); tablemapexpup[ipdf][ix] = value;};
  void SetExpDn(pdftype ipdf, int ix, double value) {tablemapexpdn[ipdf].resize(ix+1); tablemapexpdn[ipdf][ix] = value;};
  void SetModelUp(pdftype ipdf, int ix, double value) {tablemapmodelup[ipdf].resize(ix+1); tablemapmodelup[ipdf][ix] = value;};
  void SetModelDn(pdftype ipdf, int ix, double value) {tablemapmodeldn[ipdf].resize(ix+1); tablemapmodeldn[ipdf][ix] = value;};

 private:
  float Q2value;
  double Xmin;
  double Xmax;
  int    NxValues;
  int    NPdfs;
  map <pdftype, vector <double> > tablemap;
  map <pdftype, vector <double> > tablemapup;
  map <pdftype, vector <double> > tablemapdn;
  map <pdftype, vector <double> > tablemapexpup;
  map <pdftype, vector <double> > tablemapexpdn;
  map <pdftype, vector <double> > tablemapmodelup;
  map <pdftype, vector <double> > tablemapmodeldn;
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
  void pdfRotate(string dirname, string label);  //rotate PDF
  void pdfSet(string dirname, string label);  // get single set
  pdferr err;   //Type of PDF uncertainties
  bool model;   //Model PDF uncertainties
  bool par;   //Parametrisation PDF uncertainties

  map <float, Pdf> Central;               //Q2 values central PDF map with full uncertainty
//  map <float, Pdf> Exp;                   //Q2 values central PDF map with experimental uncertainty
//  map <float, Pdf> Model;                 //Q2 values central PDF map with exp + model uncertainty
  map <float, vector <Pdf> > Errors;      //Q2 values PDF errors map
  map <float, vector <Pdf> > ModelErrors; //Q2 values PDF errors map
  map <float, vector <Pdf> > ParErrors;   //Q2 values PDF errors map
  map <float, Pdf> Up;                    //Q2 values up PDF map
  map <float, Pdf> Down;                  //Q2 values down PDF map
  map <float, Pdf> UpExp;                 //Q2 values up PDF map for experimental uncertainty
  map <float, Pdf> DownExp;               //Q2 values down PDF map for experimental uncertainty
  map <float, Pdf> UpModel;               //Q2 values up PDF map for exp + model uncertainty
  map <float, Pdf> DownModel;             //Q2 values down PDF map for exp + model uncertainty

  vector <pdfshift> pdfshifts;          // PDF shifts for PDF profiling
  vector <vector <double> > cor_matrix; // correlation matrix of pdf shifts for PDF profiling
  vector<double> mcw;                   // Bayesian weights for PDF reweighting
  vector<double> mcchi2;                   // Chi2 for PDF reweighting (needed for control plots)




};

#endif
