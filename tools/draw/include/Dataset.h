#ifndef Dataset_h
#define Dataset_h

#include <vector>
#include <map>
#include <TH1F.h>

using namespace std;

//Class storing all the information of one subplot
class Subplot
{
private:
  string title, xlabel, ylabel; //data label, used for data legend
  string extralabel; //extra label
  string lumilabel;  //luminosity label
  string experiment; //experiment or collaboration
  TH1F* href;
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
  int lowrange, uprange;
  float yminr, ymaxr;
  float ymin, ymax;

  bool maketgraph;
  vector <float> valx;
  vector <float> data;
  vector <float> uncorerr;
  vector <float> toterr;
  vector <float> theory;
  vector <float> theoryshifted;
  vector <float> therrup;
  vector <float> therrdown;
  vector <float> pulls;
  vector <float> bins1;
  vector <float> bins2;

  bool hastherr;
  bool valid;
 public:
  Subplot() {};
  Subplot(string plotoptions);
  void AddPoint(map <string, float> fline);
  void Init(string label, int dataindex, int subplotindex);
  bool HasTherr() {return hastherr;};
  bool IsValid() {return valid;};
  void Draw(TH1F* histo, string opt);
  string gettitle() {return title;};
  string getextralabel() {return extralabel;};
  string getlumilabel() {return lumilabel;};
  string getexperiment() {return experiment;};
  TH1F* getref() {return href;};
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
  float getlowrange() {return lowrange;};
  float getuprange() {return uprange;};
  float getyminr() {return yminr;};
  float getymaxr() {return ymaxr;};
  float getymin() {return ymin;};
  float getymax() {return ymax;};
  bool bincenter() {return maketgraph;};
  int nbins() {return bins1.size();};
};

//Class storing all the information of one dataset
class Dataset
{
private:
  string name; //dataset name
  int index; //dataset index
  string label; //directory label, used for theory legend
public:
  Dataset() {};
  Dataset(int dataindex, string dataname);
  map <int, Subplot> subplots;
  string GetName() {return name;};
  int GetIndex() {return index;};
};


//Class storing all the datasets of one directory
class Data
{
 public:
  Data() {};
  Data(string dirname, string label);

  map <int, Dataset> datamap;      //index-dataset map
};
#endif
