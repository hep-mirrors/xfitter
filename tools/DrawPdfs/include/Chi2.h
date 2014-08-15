#ifndef Chi2_h
#define Chi2_h

#include <map>
#include <string>

using namespace std;

struct chi2type
{
  double chi2;
  double chi2_00;
  double chi2_log;
  int dof;
};

class Chi2
{
 public:
  Chi2() {};
  Chi2(string dirname, string label);
  
  //Partial, correlated and total Chi2
  map <int, chi2type> chi2list; //datasetindex-chi2 map
  chi2type chi2corr;
  chi2type chi2log;
  chi2type chi2tot;

  //Shifts (yet to be moved here)
  map <string, float> shiftlist; //shifts name-value map
};
#endif
