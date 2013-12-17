#include <vector>
#include <TH1F.h>
#include <TCanvas.h>
#include <DataSet.h>

using namespace std;

class dataseth 
{
private:
  string name; //dataset name
  string label; //directory label, used for theory legend
  string title; //data label, used for data legend
  string extralabel; //extra label
  string experiment; //experiment or collaboration
  TH1F* hdata;
  TH1F* hdatatot;
  TH1F* hth;
  TH1F* hthshift;
  TH1F* htherr;
  TH1F* htherrup;
  TH1F* htherrdown;
  TH1F* hpull;
  TH1F* r_th;
  TH1F* r_thshift;
  TH1F* r_therr;
  TH1F* r_therrup;
  TH1F* r_therrdown;
  bool logx, logy;
  float xmin, xmax;
  float yminr, ymaxr;
public:
  dataseth(string dataname, string dir, string lab, DataSet* DT, int subp);
  string getname() {return name;};
  string getlabel() {return label;};
  string gettitle() {return title;};
  string getextralabel() {return extralabel;};
  string getexperiment() {return experiment;};
  TH1F* getdata() {return hdata;};
  TH1F* getdatatot() {return hdatatot;};
  TH1F* getth() {return hth;};
  TH1F* getthshift() {return hthshift;};
  TH1F* gettherr() {return htherr;};
  TH1F* gettherrup() {return htherrup;};
  TH1F* gettherrdown() {return htherrdown;};
  TH1F* getrth() {return r_th;};
  TH1F* getrthshift() {return r_thshift;};
  TH1F* getrtherr() {return r_therr;};
  TH1F* getrtherrup() {return r_therrup;};
  TH1F* getrtherrdown() {return r_therrdown;};
  TH1F* getpull() {return hpull;};
  bool getlogx() {return logx;};
  bool getlogy() {return logy;};
  float getxmin() {return xmin;};
  float getxmax() {return xmax;};
  float getyminr() {return yminr;};
  float getymaxr() {return ymaxr;};
};


TCanvas *DataPainter(int dataindex,  vector <dataseth> datahistos);
