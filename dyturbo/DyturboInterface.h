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
  void SetOrdScales(int iord, double kmuren, double kmufac, double kmures = -1.);

  void Calculate(const double muren, const double mufac, const double mures);
  
  vector<double> upedge;
  vector<double> lowedge;
  //values of the cross section in each bin
  vector<double> values;
  //bin widths
  vector<double> binwidth;
  double yl, yh;
  double ml, mh;
  string infile;
};


#endif
