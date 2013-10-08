#ifndef CommandParser_h
#define CommandParser_h
#include <string>
#include <vector>
#include <iostream>

extern float txtsize;
extern float offset;
extern float lmarg;


using namespace std;

class CommandParser
{
 public:
  CommandParser() {};
  CommandParser(int argc, char **argv);

  //pdf options
  bool dobands, splitplots, filledbands, asymbands, logx;
  float rmin, rmax;
  float xmin, xmax;

  //data pulls options
  //  bool logx, logy;
  //  float xmin, xmax;
  bool therr, points;

  //general options
  bool pdf;
  string outdir;
  vector <string> dirs, labels;
  int colors[6];
  int styles[6];

 private:
  vector <string> allargs;
  void help(void)
  {
    cout << endl;
    cout << "program usage:" << endl;
    cout << allargs[0] << " [options] dir1[:label1] [dir2:[label2]] [...]" << endl;
    cout << endl;
    cout << "First directory is used as reference for PDF ratio plots" << endl;
    cout << "and to display data in data pulls plots." << endl;
    cout << "Directory labels are used in the legends of the plots, to add spaces and" << endl;
    cout << "special characters to the labels use quotation marks '' (ex. dir:'HERA I')." << endl;
    cout << "To specify greek letters and latex commands in the labels use the ROOT notation" << endl;
    cout << "(#alpha #bar{u})." << endl;
    cout << "It is possible to specify up to 6 directories, you need to specify at least one directory." << endl;
    cout << endl;
    cout << "general options:" << endl;
    cout << "\t --help" << endl;
    cout << "\t \t Show this help" << endl;
    cout << "\t --pdf (require ps2pdf)" << endl;
    cout << "\t \t Produce and addition plot file in pdf format" << endl;
    cout << "\t --colorpattern <1-3>" << endl;
    cout << "\t \t Select among 3 additional color patterns" << endl;
    cout << "options for pdf plots:" << endl;
    cout << "\t --bands" << endl;
    cout << "\t \t Draw PDF uncertainty band" << endl;
    cout << "\t --asymbands" << endl;
    cout << "\t \t PDF bands are not symmetrised" << endl;
    cout << "\t --outdir <output directory>" << endl;
    cout << "\t \t Specify output directory" << endl;
    cout << "\t --splitplots" << endl;
    cout << "\t \t Produce also additional eps files for each plot" << endl;
    cout << "\t --filledbands" << endl;
    cout << "\t \t Filled uncertainty bands, usefull for sensitivity studies" << endl;
    cout << "\t --ratiorange min:max" << endl;
    cout << "\t \t Specify y axis range in PDF ratio plots" << endl;
    cout << "\t --xrange min:max" << endl;
    cout << "\t \t Specify x axis range in PDF plots: minimum x is 0.0001, maximum x is 1" << endl;
    cout << "\t --no-logx" << endl;
    cout << "\t \t Linear x scale in PDF plots" << endl;
    cout << "options for data pulls plots:" << endl;
    cout << "\t --therr" << endl;
    cout << "\t \t Plot theory errors if availables" << endl;
    cout << "\t --points" << endl;
    cout << "\t \t Plot theory as tilted marker points (with vertical error bars) instead of continous lines (with dashed error area)" << endl;
    cout << "\t to set axis titles, axis range and log scales add PlotDesc options in the data file" << endl;
    cout << "\t Example:" << endl;
    cout << "\t &PlotDesc" << endl;
    cout << "\t    PlotN = 1" << endl;
    cout << "\t    PlotDefColumn = 'y2'" << endl;
    cout << "\t    PlotDefValue = 0., 5." << endl;
    cout << "\t    PlotOptions(1)  = 'XTitle: y_{Z} [GeV/c] @YTitle: d#sigma/dy_{Z} [pb]  @Xmin:0.0@Xmax:2.5@Ylog'" << endl;
    cout << "\t &End" << endl;
    cout << endl;
  };  
};

extern CommandParser opts;
#endif
