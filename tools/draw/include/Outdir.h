#ifndef  Outdir_h
#define  Outdir_h

#include "PdfData.h"
#include "Dataset.h"
#include "Chi2.h"
#include "Par.h"
#include "Chi2scanData.h"

#include <TGraphAsymmErrors.h>

#include <vector>
#include <map>
#include <string>

extern map <string, PdfData> pdfmap;  //dir label-pdf map	   
//extern map <string, PdfData> pdfmapProfiled;  //dir label-pdf map for profiled PDFs too
extern map <string, Data> datamap;    //dir label-dataset map
extern map <string, Chi2> chi2map;    //dir label-chi2 map
extern map <string, Par> parmap;      //dir label-parameters map
extern map <string, Chi2scanData> chi2scanmap; //dir label-chi2scan map

extern vector <TGraphAsymmErrors*> allgraphs; //global TGraph vector to store graphs in plots.root

using namespace std;

class  Outdir {
 private:    
  string label;
  string dirname;
  bool MCreplica, cl68, cl90, median, asym, bands3;
  bool profiled;
  bool reweighted;
  bool GK;
  bool BAY;
  bool scale68;
  bool rotated;
  int RotatedSet;
  bool singlePDF;
  int singlePDFlabel;
 public:
  Outdir() {};
  Outdir(string dir);
  string GetName() {return dirname;};
  string GetLabel() {return label;};
  bool IsMCreplica() {return MCreplica;};
  bool IsProfiled() {return profiled;};
  bool IsReweighted() {return reweighted;};
  bool IsGK() {return GK;};
  bool IsBAY() {return BAY;};
  bool Scale68() {return scale68;};
  bool IsRotated() {return rotated;};
  bool IsSingleSet() {return singlePDF; };
  int  rSet()  {return RotatedSet; };
  int  pdfSet() {return singlePDFlabel;}; 
  bool Is68cl() {return cl68;};
  bool Is90cl() {return cl90;};
  bool IsMedian() {return median;};
  bool IsAsym() {return asym;};
  bool Is3bands() {return bands3;};
  vector <string> dirlist;
};

extern map <string, Outdir> outdirs; //label-directory map

extern vector <float> q2list();
extern vector <int> datalist();
extern vector <int> chi2list();
extern vector <int> parlist();
extern vector <string> chi2scanlist();
extern int finddataindex(string name);
extern string finddataname(int index);
extern string findparname(int index);

#endif
