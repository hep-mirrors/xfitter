#ifndef Par_h
#define Par_h

#include <string>
#include <map>
#include <vector>

using namespace std;

extern string fitstat(string dir);
extern string hessestat(string dir);

struct partype
{
  string name;
  double value;
  double error_p;
  double error_m;
};

class Par
{
 public:
  Par() {};
  Par(string dirname, string label);
  
  //Parameters map
  map <int, partype> parlist; //parindex-partype map

  //Parameters vector to store MC replica
  map <int, vector <double> > MCparams;
  string fitstatus;
  string uncertainties;
};

#endif
