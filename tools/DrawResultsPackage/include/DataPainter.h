#include <vector>
#include <TH1F.h>
#include <TCanvas.h>

using namespace std;

class dataseth 
{
private:
  string name; //dataset name
  string label; //directory label
  TH1F* hdata;
  TH1F* hdatatot;
  TH1F* hth;
  TH1F* hthshift;
  TH1F* hpull;
public:
  dataseth(string dataname, string dir, string lab,
	   vector <float> bins1, vector <float> bins2, 
	   vector <float> data, vector <float> uncorerr, vector <float> toterr, 
	   vector <float> theory, vector <float> theoryshifted, 
	   vector <float> pulls);

  string getname() {return name;};
  string getlabel() {return label;};
  TH1F* getdata() {return hdata;};
  TH1F* getdatatot() {return hdatatot;};
  TH1F* getth() {return hth;};
  TH1F* getthshift() {return hthshift;};
  TH1F* getpull() {return hpull;};
};


TCanvas *DataPainter(int dataindex,  vector <dataseth> datahistos);
