#ifndef  Outdir_h
#define  Outdir_h

#include <vector>
#include <map>
#include <string>

#include "PdfData.h"
#include "Dataset.h"
#include "Chi2.h"
#include "Par.h"


extern map <string, PdfData> pdfmap;  //dir label-pdf map	   
extern map <string, Data> datamap;    //dir label-dataset map
extern map <string, Chi2> chi2map;    //dir label-chi2 map
extern map <string, Par> parmap;      //dir label-parameters map

using namespace std;

class  Outdir {
 private:    
  string label;
  string dirname;
  bool MCreplica, cl68, cl90, median, asym;

 public:
  Outdir() {};
  Outdir(string dir);
  string GetName() {return dirname;};
  string GetLabel() {return label;};
  bool IsMCreplica() {return MCreplica;};
  bool Is68cl() {return cl68;};
  bool Is90cl() {return cl90;};
  bool IsMedian() {return median;};
  bool IsAsym() {return asym;};
  vector <string> dirlist;
};

extern map <string, Outdir> outdirs; //label-directory map

extern vector <float> q2list();
extern vector <int> datalist();
extern vector <int> chi2list();
extern vector <int> parlist();
extern int finddataindex(string name);
extern string finddataname(int index);
extern string findparname(int index);

#endif
