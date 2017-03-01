#ifndef CommandParser_h
#define CommandParser_h
#include <string>
#include <vector>
#include <map>
#include <iostream>

extern float txtsize;
extern float offset;
extern float lmarg;
extern float rmarg;
extern float tmarg;
extern float bmarg;
extern float marg0;

using namespace std;

//Service functions
extern vector<string> Round(double value, double error = 0, bool sign = false);

class CommandParser
{
 public:
  CommandParser() {};
  CommandParser(int argc, char **argv);

  //pdf options
  bool dobands, filledbands, asym, logx;
  float rmin, rmax;
  double xmin, xmax;
  bool abserror, relerror;
  bool q2all;
  bool cl68, cl90, median;
  int plotsperpage;
  bool scale68;

  //data pulls options
  bool therr, points;
  string theorylabel;
  bool twopanels,threepanels;
  bool multitheory;
  bool nothshifts;
  bool onlytheory;
  bool threlerr;
  bool ratiototheory;
  bool diff;
  bool noupband;
  int errbandcol;
  
  //shifts options
  int spp, shgth;
  bool adjshift;

  //tables options
  bool chi2nopdf;
  bool logpenalty;
  string font;

  //general options
  bool splitplots;
  string ext;
  string format;
  bool root;
  string outdir;
  vector <string> dirs, labels;
  map <string, int> colors;
  map <string, int> styles;
  map <string, int> lstyles;
  map <string, int> markers;
  int col[6];
  int styl[6];
  int lstyl[6];
  int mark[6];
  float lwidth;
  float resolution, pagewidth;
  bool nodata;
  bool nopdfs;
  bool noshifts;
  bool notables;
  bool version, drawlogo;
  bool cms, cmspreliminary;
  bool atlas, atlaspreliminary, atlasinternal;
  bool cdfiipreliminary;
  bool profile, reweight, BAYweight, GKweight;
  bool bw;
  bool looseRepSelection;
    
private:
  vector <string> allargs;
  void help(void)
  {
    cout << endl;
    cout << "program usage:" << endl;
    cout << allargs[0] << " [options] dir1[:label1] [dir2:[label2]] [MC:dirpattern:[label3]] [...]" << endl;
    cout << endl;
    cout << "First directory is used as reference for PDF ratio plots" << endl;
    cout << "and to display data in data plots." << endl;
    cout << "Directory labels are used in the legends of the plots, to add spaces and" << endl;
    cout << "special characters to the labels use quotation marks '' (ex. dir:'HERA I')." << endl;
    cout << "To specify greek letters and latex commands in the labels use the ROOT notation" << endl;
    cout << "(#alpha #bar{u})." << endl;
    cout << "It is possible to specify up to 6 directories, you need to specify at least one directory." << endl;
    cout << endl;
    cout << "Monte Carlo replica directories:" << endl;
    cout << "\t To specify a pattern of directories containing Monte Carlo replica" << endl;
    cout << "\t use the prefix \"MC:\" as in \"MC:dirpattern\"." << endl;
    cout << "\t If \"dirpattern\" is a directory all the subdirectories" << endl;
    cout << "\t are considered as Monte Carlo replica runs." << endl;
    cout << "\t If \"dirpattern\" is not a directory" << endl;
    cout << "\t all the directories \"dirpattern*\" are considered as Monte Carlo replica runs" << endl;
    cout << "general options:" << endl;
    cout << "\t --help" << endl;
    cout << "\t \t Show this help" << endl;
    cout << "\t --outdir <output directory>" << endl;
    cout << "\t \t Specify output directory" << endl;
    cout << "\t --eps" << endl;
    cout << "\t \t Plots are saved in eps format" << endl;
    cout << "\t --root" << endl;
    cout << "\t \t Save all the plots in a root file" << endl;
    cout << "\t --splitplots-eps" << endl;
    cout << "\t \t Produce also additional eps files for each plot" << endl;
    cout << "\t --splitplots-pdf" << endl;
    cout << "\t \t Produce also additional eps and pdf files for each plot" << endl;
    cout << "\t --splitplots-png" << endl;
    cout << "\t \t Produce additional png files for each plot" << endl;
    cout << "\t --colorpattern <1-3>" << endl;
    cout << "\t \t Select among 3 additional color patterns" << endl;
    cout << "\t --thicklines" << endl;
    cout << "\t \t Thicker lines in all plots (better for small plots in slides)" << endl;
    cout << "\t --largetext" << endl;
    cout << "\t \t Larger text size in all plots (better for small plots in slides)" << endl;
    cout << "\t --bw" << endl;
    cout << "\t \t Black and white compatible plots" << endl;
    cout << "\t --lowres" << endl;
    cout << "\t \t Low resolution logo (smaller file)" << endl;
    cout << "\t --highres" << endl;
    cout << "\t \t High resolution logo (paper quality)" << endl;
    cout << "\t --no-version" << endl;
    cout << "\t \t Do not show version on logo" << endl;
    cout << "\t --plots-per-page <N>" << endl;
    cout << "\t \t Number of rows and columns of PDF and data plots per page, default value is 2" << endl;
    cout << "\t --loose-mc-replica-selection" <<endl;
    cout << "\t \t Do not check for fit convergence for MC replica " <<endl;
    

    cout << "options for pdf plots:" << endl;
    cout << "\t --no-pdfs" << endl;
    cout << "\t \t PDF plots are not produced" << endl;
    cout << "\t --bands" << endl;
    cout << "\t \t Draw PDF uncertainty bands" << endl;
    cout << "\t --profile" << endl;
    cout << "\t \t Draw Profiled PDF (only for Hessian sets)" << endl;
    cout << "\t \t To set this option only for one directory, use the syntax profiled:directory[:label]" << endl;
    cout << "\t Example: xfitter-draw profile:output:\"profiled\" output:\"not-profiled\"" << endl;
    cout << "\t --reweight(-BAY/-GK)" << endl;
    cout << "\t \t Draw Reweighted PDF (only for MC replica sets)" << endl;
    cout << "\t \t To set this option only for one directory, use the syntax reweight-BAY:directory[:label]" << endl;
    cout << "\t \t To use the Giele-Keller weights instead of Bayesian weights, use the syntax reweight-GK:directory[:label]" << endl;
    cout << "\t Example: xfitter-draw reweight-BAY:output:\"reweighted\" output:\"not-reweighted\"" << endl;

    cout << "\t options for rotation:" << endl;
    cout << "\t \t Draw Rotated PDF (only for Hessian sets)" << endl;
    cout << "\t \t To set this option, use the syntax rotate:<n>:directory[:label]" << endl;
    cout << "\t Example: xfitter-draw rotate:5:output:\"rotated-5\" output:\"not-rotated\"" << endl;

    cout << "\t options for individual eigen sets:" << endl;
    cout << "\t \t Draw individual eigen (or MC replica) set from a complete run" << endl;
    cout << "\t \t To set this option, use the syntax set:<n>:directory[:label]" << endl;
    cout << "\t Example: xfitter-draw --bands set:5:output:\"set 5\" output:\"Total Band\"" << endl;


    cout << "\t --filledbands" << endl;
    cout << "\t \t Filled uncertainty bands, usefull for sensitivity studies" << endl;
    cout << "\t --ratiorange min:max" << endl;
    cout << "\t \t Specify y axis range in PDF ratio plots" << endl;
    cout << "\t --xrange min:max" << endl;
    cout << "\t \t Specify x axis range in PDF plots," << endl;
    cout << "\t \t default minimum and maximum x are determined by settings" << endl;
    cout << "\t \t in the first reference directory" << endl;
    cout << "\t --no-logx" << endl;
    cout << "\t \t Linear x scale in PDF plots" << endl;
    cout << "\t --absolute-errors" << endl;
    cout << "\t \t Plot absolute pdf uncertainties centered around 0 in PDF ratio plots" << endl;
    cout << "\t --relative-errors" << endl;
    cout << "\t \t Plot relative pdf uncertainties centered around 1 in PDF ratio plots" << endl;
    cout << "\t --q2all" << endl;
    cout << "\t \t Plot PDF at all stored values of Q2. By default PDF are plotted only at the starting scale Q0" << endl;
    cout << "options for data plots:" << endl;
    cout << "\t --no-data" << endl;
    cout << "\t \t Data plots are not produced" << endl;
    cout << "\t --therr" << endl;
    cout << "\t \t Plot theory uncertainties if available" << endl;
    cout << "\t --noupband" << endl;
    cout << "\t \t Do not plot theory uncertainties in the upper panel" << endl;
    cout << "\t --nothshifts" << endl;
    cout << "\t \t Do not plot theory+shifts lines" << endl;
    cout << "\t --points" << endl;
    cout << "\t \t Plot theory as displaced marker points (with vertical error bars) instead of continous lines (with dashed error area)" << endl;
    cout << "\t --theory <label>" << endl;
    cout << "\t \t Change the \"Theory\" legend to the specified label" << endl;
    cout << "\t --2panels" << endl;
    cout << "\t \t Plot additional right bottom panels with pulls" << endl;
    cout << "\t --3panels" << endl;
    cout << "\t \t Plot additional right mid panels with theory+shifts" << endl;
    cout << "\t --only-theory" << endl;
    cout << "\t \t Do not plot data" << endl;
    cout << "\t --ratio-to-theory" << endl;
    cout << "\t \t Use theory as reference for ratio plots" << endl;
    cout << "\t --theory-rel-errors" << endl;
    cout << "\t \t Do not plot data, use theory as reference for ratio plots, and plot relative theory uncertainties" << endl;
    cout << "\t --diff" << endl;
    cout << "\t \t Plot difference of theory-data instead of ratio theory/data" << endl;
    cout << "\t --greenband" << endl;
    cout << "\t \t The total experimental uncertainty is shown with a green band" << endl;
    cout << "\t --blueband" << endl;
    cout << "\t \t The total experimental uncertainty is shown with a blue band" << endl;
    cout << "\t --multitheory" << endl;
    cout << "\t \t Plot ratios of theory predictions in separate panels (up to three)" << endl;
    cout << "options for shifts plots:" << endl;
    cout << "\t --no-shifts" << endl;
    cout << "\t \t Shifts plots are not produced" << endl;
    cout << "\t --shifts-per-plot <N>" << endl;
    cout << "\t \t Number of shifts shown in each plot, default is 30, maximum is 40" << endl;
    cout << "\t --shifts-heigth <N>" << endl;
    cout << "\t \t Heigth reserved for each shift in points, minimum is 20, maximum is 200" << endl;
    cout << "options for tables:" << endl;
    cout << "\t --no-tables" << endl;
    cout << "\t \t Chi2 and parameter tables are not produced" << endl;
    cout << "\t --chi2-nopdf-uncertainties" << endl;
    cout << "\t \t When chi2 is evaluated with the LHAPDFError routine, this option adds to " << endl;
    cout << "\t \t the chi2 table the chi2 evaluated without PDF uncertainties separated by a | character" << endl;
    cout << "\t --partial-log-penalty" << endl;
    cout << "\t \t Show within brackets the partial log penalty chi2 for each dataset in the chi2 table" << endl;
    cout << "\t \t This option and the above --chi2-nopdf-uncertainties are mutually exclusive" << endl;
    cout << "\t --helvet-fonts" << endl;
    cout << "\t \t Use helvetica fonts in tables (default is palatino)" << endl;
    cout << "\t --cmbright-fonts" << endl;
    cout << "\t \t Use Computer Modern Bright fonts in tables (default is palatino)" << endl;
    cout << "Statistical option for PDF error bands in PDF plots and parameter errors in parameter table." << endl;
    cout << "\t The following options apply to MC-replica and MC error PDF." << endl;
    cout << "\t The \"asym\" option applies also to asymmetric hessian error PDF." << endl;
    cout << "\t All the options can be set globally, or separately for each directory" << endl;
    cout << "\t (or for each pattern of directories in the case of MC replica runs)." << endl;
    cout << "\t To set the options for a directory, use the syntax [MC:]<option1:>[option2:]directory[:label]" << endl;
    cout << "\t Example: xfitter-draw MC:68cl:asym:MYMCReplicaRuns" << endl;
    cout << "\t Directory specific options take precedence with respect to global options." << endl;
    cout << "\t --median" << endl;
    cout << "\t \t Use median instead of average for the central values of PDF and parameters" << endl;
    cout << "\t --68cl" << endl;
    cout << "\t \t Evaluate 68 cl for PDF and parameters uncertainties. The option --median is activated automatically" << endl;
    cout << "\t \t Can be use in conjunction with --asymbands to get asymmetric errors" << endl;
    cout << "\t --90cl" << endl;
    cout << "\t \t Evaluate 90 cl for PDF and parameters uncertainties. The option --median is activated automatically" << endl;
    cout << "\t --asym" << endl;
    cout << "\t \t Evaluate asymmetric errors when possible" << endl;
    cout << "\t \t Can be use in conjunction with --asymbands to get asymmetric errors" << endl;
    cout << "\t --scale68" << endl;
    cout << "\t \t Scale PDF uncertainty bands by 1.645 (affects PDF and data plots)" << endl;
    cout << endl;
    cout << "\t to set axis titles, axis range and log scales add PlotDesc options in the data file" << endl;
    cout << "\t Example:" << endl;
    cout << "\t &PlotDesc" << endl;
    cout << "\t    PlotN = 1 ! SubPlot index" << endl;
    cout << "\t    PlotDefColumn = 'y2' ! Variable used to divide the data in SubPlots" << endl;
    cout << "\t    PlotDefValue = 0., 5.! Ranges of PlotDefColumn used to divide the data in SubPlots" << endl;
    cout << "\t    PlotVarColumn = 'x'! Variable providing bin center information (use only if bin edges are missing)" << endl;
    cout << "\t    PlotOptions(1)  = 'Experiment:ATLAS@Title: pp #rightarrow Z@XTitle: y_{Z} @YTitle: d#sigma/dy_{Z} [pb]  @ExtraLabel:#int L = 100 fb^{-1}@Xmin:0.0@Xmax:2.5@Xlog@Ylog@YminR:0.91@YmaxR:1.09'" 
	 << endl;
    cout << "\t &End" << endl;
    cout << endl;
    cout << "Option for 3-band PDF uncertainty bands (HERAPDF style) in PDF plots." << endl;
    cout << "\t --bands 3bands:<dir-full-uncert>" << endl;
    cout << "\t \t draw PDFs with three uncertainty bands: experimental (red), model (yellow) parametrisation (green)." << endl;
    cout << "\t The model uncertainty originates from variation of model parameters (e.g. masses of heavy quarks, Q^2 cuts on data, etc.),"  << endl; 
    cout << "\t parametrisation - from variations of the parameters in the fit and variation of the starting scale Q_0^2."<< endl;
    cout << " \t Directory <dir-full-uncert> should have fit results for experimental, model and parametrisation variations. " << endl;
    cout << " \t The file names for experimental variations should follow convention as follows: " << endl;
    cout << " \t pdfs_q2val_s01m_01.txt  " << endl;
    cout << " \t pdfs_q2val_s01p_01.txt  " << endl;
    cout << " \t ...  " << endl;
    cout << " \t where s01m stands for experimental error (s), number of fitted parameters (01 to N) and minus (m) or plus (p) variation; " << endl;
    cout << " \t the last number stands for the index of the Q^2 value at which PDFs are drawn (defined in Q2VAL in steering.txt). " << endl;
    cout << " \t Similarly, m11m stands for model uncertainty and the number should start from N+1 (here assuming that N=10 for exp errors)." << endl;
    cout << " \t Finally, p14s stands for parametrisation uncertainty and the number should start from N+K+1 (here assuming that K=3 for model errors)." << endl;
    cout << " \t NOTE: if command '--bands <dir-full-uncert>' is used, the total uncertainty in red is drawn." << endl;
    cout << endl;
  };  
};

extern CommandParser opts;
#endif
