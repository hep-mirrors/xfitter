#ifndef Chi2scanData_h
#define Chi2scanData_h

#include <string>
#include <map>

using namespace std;

class Chi2scanData
{
 public:
  Chi2scanData() {};
  Chi2scanData(string dirname);
  string label;
  double min, delta, chi2min;
  map <double, double> chi2;

  map <int, double> pdfmin;
  map <int, double> pdfdelta;
  map <int, double> pdfchi2min;
  map <int, map <double, double> > pdfchi2;
};
#endif
