#ifndef Par_h
#define Par_h

#include <string>
#include <map>

using namespace std;
struct partype
{
  string name;
  double value;
  double error;
};

class Par
{
 public:
  Par() {};
  Par(string dirname, string label);
  
  //Parameters map
  map <int, partype> parlist; //parindex-partype map
  string fitstatus;
};

#endif
