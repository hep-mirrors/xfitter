#ifndef DyturboInterface_h
#define DyturboInterface_h

#include <string>
#include <vector>

using namespace std;

class Dyturbo
{
 public:
  Dyturbo(string file);

  void SetBins(vector <double> ledge, vector <double> uedge, double ylow, double yhigh, double mlow, double mhigh);
  void SetOrdScales(int iord, double kmuren, double kmufac, double kmures = -1., double muC3 = 1.);

  void Calculate(const double muren, const double mufac, const double mures, const double muC3 = 1.);
  
  vector<double> upedge;
  vector<double> lowedge;
  //values of the cross section in each bin
  vector<double> values;
  //bin widths
  vector<double> binwidth;
  double yl, yh;
  double ml, mh;
  int proc;
  string infile;
  int order;

  static string pdfname;
  static int pdfmember;
};


#endif
