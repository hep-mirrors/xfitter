#ifndef CommandParser_h
#define CommandParser_h
#include <string>
#include <vector>

using namespace std;

class CommandParser
{
 public:
  CommandParser() {};
  CommandParser(int argc, char **argv);
  bool dobands, splitplots, filledbands;
  string outdir;
  
  vector <string> dirs, labels;

  float rmin, rmax;

  int colors[6];
  int styles[6];

 private:
  
};

extern CommandParser opts;
#endif
