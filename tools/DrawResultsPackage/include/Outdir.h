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
 protected:    
  string label;
  string dirname;
  
 public:
  Outdir(string dir, string lab);
  string GetName() {return dirname;}
  string GetLabel() {return label;}
};

extern vector <float> q2list();
extern vector <int> datalist();
extern vector <int> chi2list();
extern vector <int> parlist();
extern int finddataindex(string name);
extern string finddataname(int index);
extern string findparname(int index);

#endif
