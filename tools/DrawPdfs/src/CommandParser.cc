#include "CommandParser.h"

#include <TStyle.h>
#include <TMath.h>

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <algorithm>

float txtsize = 0.043;
float offset = 1.5;
float lmarg = 0.15;
float rmarg = 0.05;
float tmarg = 0.05;
float bmarg = 0.1;
float marg0 = 0.003;

CommandParser::CommandParser(int argc, char **argv):
  dobands(false),
  asym(false),
  logx(true),
  filledbands(false),
  rmin(0),
  rmax(0),
  xmin(-1),
  xmax(-1),
  relerror(false),
  abserror(false),
  q2all(false),
  cl68(false),
  cl90(false),
  median(false),
  splitplots(false),
  root(false),
  format("pdf"),
  ext("eps"),
  resolution(800),
  pagewidth(20),
  lwidth(2),
  therr(false),
  points(false),
  theorylabel("Theory"),
  onlytheory(false),
  ratiototheory(false),
  diff(false),
  twopanels(false),
  threepanels(false),
  version(true),
  drawlogo(true),
  nodata(false),
  nopdfs(false),
  noshifts(false),
  notables(false),
  spp(30),
  shgth(40),
  adjshift(true),
  chi2nopdf(false),
  font("palatino"), //font("helvet"), //font("modernbright"),
  cms(false),
  cmspreliminary(false),
  atlas(false),
  atlaspreliminary(false),
  atlasinternal(false),
  cdfiipreliminary(false),
  outdir("")
{
  //initialise colors and styles
  col[0] = kRed + 2;
  col[1] = kBlue + 2;
  col[2] = kGreen + 3;
  col[3] = kOrange + 7;
  col[4] = kAzure + 1;
  col[5] = kMagenta + 1;

  styl[0] = 3354;
  styl[1] = 3345;
  styl[2] = 3359;
  styl[3] = 3350;
  styl[4] = 3016;
  styl[5] = 3020;

  mark[0] = 24;
  mark[1] = 25;
  mark[2] = 26;
  mark[3] = 32;
  mark[4] = 31;
  mark[5] = 27;

  //Set Hatches style
  gStyle->SetHatchesSpacing(2);
  gStyle->SetHatchesLineWidth(lwidth);

  //read all command line arguments
  for (int iar = 0; iar < argc; iar++)
    allargs.push_back(argv[iar]);

  //search for options
  for (vector<string>::iterator it = allargs.begin() + 1; it != allargs.end(); it++)
    if ((*it).find("--") == 0)
      {
	if (*it == "--help")
	  {
	    help();
	    exit(0);
	  }
	else if (*it == "--thicklines")
	  lwidth = 3;
	else if (*it == "--lowres")
	  {
	    resolution = 400;
	    //	    pagewidth = 10;
	  }
	else if (*it == "--highres")
	  {
	    resolution = 2400;
	    //	    pagewidth = 60;
	  }
	else if (*it == "--no-version")
	  version = false;
	else if (*it == "--no-logo")
	  drawlogo = false;
	else if (*it == "--no-data")
	  nodata = true;
	else if (*it == "--no-pdfs")
	  nopdfs = true;
	else if (*it == "--no-shifts")
	  noshifts = true;
	else if (*it == "--no-tables")
	  notables = true;
	else if (*it == "--chi2-nopdf-uncertainties")
	  chi2nopdf = true;
	else if (*it == "--helvet-fonts")
	  font = "helvet";
	else if (*it == "--cmbright-fonts")
	  font = "modernbright";
	else if (*it == "--shifts-per-page")
	  {
	    adjshift = false;
	    spp = atoi((*(it+1)).c_str());
	    spp = max(1, spp);
	    spp = min(40, spp);
	    allargs.erase(it+1);
	  }
	else if (*it == "--shifts-heigth")
	  {
	    adjshift = false;
	    shgth = atoi((*(it+1)).c_str());
	    shgth = max(20, shgth);
	    shgth = min(200, shgth);
	    allargs.erase(it+1);
	  }
	else if (*it == "--cms")
	  cms = true;
	else if (*it == "--cms-preliminary")
	  cmspreliminary = true;
	else if (*it == "--atlas")
	  atlas = true;
	else if (*it == "--atlas-internal")
	  atlasinternal = true;
	else if (*it == "--atlas-preliminary")
	  atlaspreliminary = true;
	else if (*it == "--cdfii-preliminary")
	  cdfiipreliminary = true;
	else if (*it == "--hidden")
	  {
	    cout << endl;
	    cout << "Hidden options" << endl;
	    cout << "Please use this options only if you are authorised from your collaboration to do so" << endl;
	    cout <<  "--cms" << endl;
	    cout <<  "--cms-preliminary" << endl;
	    cout <<  "--atlas" << endl;
	    cout <<  "--atlas-internal" << endl;
	    cout <<  "--atlas-preliminary" << endl;
	    cout <<  "--cdfii-preliminary" << endl;
	    cout <<  "--no-logo" << endl;
	    cout << endl;
	    exit(-1);
	  }
	else if (*it == "--bands")
	  dobands = true;
	else if (*it == "--asym")
	  {
	    dobands = true;
	    asym = true;
	  }
	else if (*it == "--median")
	  median = true;
	else if (*it == "--68cl")
	  {
	    if (cl90 == true)
	      {
		cout << "Options --68cl and --90cl are mutually exclusive, cannot use both" << endl;
		exit(1);
	      }
	    cl68 = true;
	    median = true;
	  }
	else if (*it == "--90cl")
	  {
	    if (cl68 == true)
	      {
		cout << "Options --68cl and --90cl are mutually exclusive, cannot use both" << endl;
		exit(1);
	      }
	    cl90 = true;
	    median = true;
	  }
	else if (*it == "--absolute-errors")
	  {
	    dobands = true;
	    abserror = true;
	  }
	else if (*it == "--relative-errors")
	  {
	    dobands = true;
	    relerror = true;
	  }
	else if (*it == "--no-logx")
	  logx = false;
	else if (*it == "--q2all")
	  q2all = true;
	else if (*it == "--outdir")
	  {
	    outdir = *(it+1);
	    allargs.erase(it+1);
	  }
	else if (*it == "--eps")
	  format = "eps";
	else if (*it == "--root")
	  root = true;
	else if (*it == "--splitplots-eps")
	  {
	    splitplots = true;
	    ext = "eps";
	  }
	else if (*it == "--splitplots-pdf")
	  {
	    splitplots = true;
	    ext = "pdf";
	  }
	else if (*it == "--splitplots-png")
	  {
	    splitplots = true;
	    ext = "png";
	  }
	else if (*it == "--filledbands")
	  filledbands = true;
	else if (*it == "--ratiorange")
	  {
	    rmin = atof((*(it+1)).substr(0, (*(it+1)).find(":")).c_str());
	    rmax = atof((*(it+1)).substr((*(it+1)).find(":") + 1, (*(it+1)).size() - (*(it+1)).find(":") - 1).c_str());
	    allargs.erase(it+1);
	  }
	else if (*it == "--xrange")
	  {
	    xmin = max(0.000000000000001, atof((*(it+1)).substr(0, (*(it+1)).find(":")).c_str()));
	    xmax = min(1., atof((*(it+1)).substr((*(it+1)).find(":") + 1, (*(it+1)).size() - (*(it+1)).find(":") - 1).c_str()));
	    allargs.erase(it+1);
	  }
	else if (*it == "--colorpattern")
	  {
	    int pattern = atoi((*(it+1)).c_str());
	    if (pattern == 1)
	      {
		col[0] = kBlue + 2;
		col[1] = kYellow - 7;
		col[2] = kGreen - 3;
		col[3] = kRed + 1;
		col[5] = kOrange + 7;
		col[5] = kCyan + 1;
	      }
	    else if (pattern == 2)
	      {
		col[0] = kBlue + 2;
		col[1] = kRed + 1;
		col[2] = kYellow - 7;
		col[3] = kOrange + 7;
		col[4] = kMagenta + 1;
		col[5] = kCyan + 1;
	      }
	    else if (pattern == 3)
	      {
		col[0] = kBlue + 2;
		col[1] = kMagenta + 1;
		col[2] = kCyan + 1;
		col[3] = kRed + 1;
		col[4] = kGreen + 2;
		col[5] = kYellow + 1;
	      }

	      allargs.erase(it+1);
	  }
	else if (*it == "--therr")
	  therr = true;
	else if (*it == "--points")
	  points = true;
	else if (*it == "--theory")
	  {
	    theorylabel = *(it+1);
	    allargs.erase(it+1);
	  }
	else if (*it == "--only-theory")
	  {
	    onlytheory = true;
	    ratiototheory = true;
	  }
	else if (*it == "--ratio-to-theory")
	  ratiototheory = true;
	else if (*it == "--diff")
	  diff = true;
	else if (*it == "--2panels")
	  twopanels = true;
	else if (*it == "--3panels")
	  threepanels = true;
	else
	  {
	    cout << endl;
	    cout << "Invalid option " << *it << endl;
	    cout << allargs[0] << " --help for help " << endl;
	    cout << endl;
	    exit(-1);
	  }
	allargs.erase(it);
	it = allargs.begin();
      }
  
  for (vector<string>::iterator it = allargs.begin() + 1; it != allargs.end(); it++)
    dirs.push_back(*it);

  if (dirs.size() == 0)
    {
      cout << endl;
      cout << "Please specify at least one directory" << endl;
      cout << allargs[0] << " --help for help " << endl;
      cout << endl;
      exit(-1);
    }

  if (dirs.size() > 6)
    {
      cout << endl;
      cout << "Maximum number of directories is 6" << endl;
      cout << allargs[0] << " --help for help " << endl;
      cout << endl;
      exit(-1);
    }
}

CommandParser opts;

//Service functions
vector<string> Round(double value, double error)
{
  vector <string> result;

  int decimal = 0;

  //If no error, value is rounded to two significant digits
  if (error == 0)
    error = value;

  if (error != 0)
    decimal = -log10(fabs(error)) + 2;
  decimal = max(0, decimal);

  char Dec[2];
  sprintf (Dec, "%d", decimal);
  string D = Dec;

  char Numb[50];
  sprintf (Numb, ((string)"%." + D + "f").c_str(), value);
  result.push_back(Numb);
  sprintf (Numb, ((string)"%." + D + "f").c_str(), error);
  result.push_back(Numb);
  
  return result;
}

double Median(vector <double> xi)
{
  int n = xi.size();
  double x[n];
  std::copy (xi.begin(), xi.end()-1, x);

  return TMath::Median(n, x);
}

double cl(int sigma)
{
  switch (sigma)
    {
    case 1:
      return 0.682689492137086;
      break;
    case 2:
      return 0.954499736103642;
      break;
    case 3:
      return 0.997300203936740;
      break;
    default:
      cout << "Confidence Level interval available only for sigma = 1, 2, 3, requested: " << sigma << " sigma" << endl;
      exit(1);
    }
}

double delta(vector <double> xi, double central, double ConfLevel)
{
  double delta = 0;

  vector <double> deltaxi;
  vector<double>::iterator i = xi.begin();
  while (i != xi.end())
    {
      deltaxi.push_back(fabs(*i - central));
      ++i;
    }

  sort(deltaxi.begin(), deltaxi.end());

  vector<double>::iterator di = deltaxi.begin();
  while (di != deltaxi.end())
    {
      delta = *di;
      int index = di - deltaxi.begin() + 1;
      double prob = (double)index / (double)deltaxi.size();
      //      cout << index << "  " << *di << "  " << prob << endl;
      if (prob > ConfLevel)
	break;
      ++di;
    }
  return delta;
}
double deltaasym(vector <double> xi, double central, double& delta_p, double& delta_m, double ConfLevel)
{
  delta_m = delta_p = 0;

  vector <double> deltaxi;
  vector<double>::iterator i = xi.begin();
  while (i != xi.end())
    {
      deltaxi.push_back(*i - central);
      ++i;
    }

  sort(deltaxi.begin(), deltaxi.end());

  vector<double>::iterator di = deltaxi.begin();
  while (di != deltaxi.end())
    {
      delta_m = fabs(*di);
      int index = di - deltaxi.begin() + 1;
      double prob = (double)index / (double)deltaxi.size();
      //      cout << index << "  " << *di << "  " << prob << endl;
      if (prob >= ((1-ConfLevel)/2.))
	break;
      ++di;
    }

  di = deltaxi.end();
  while (di != deltaxi.begin())
    {
      delta_p = *di;
      int index = di - deltaxi.begin() + 1;
      double prob = (double)index / (double)deltaxi.size();
      //      cout << index << "  " << *di << "  " << prob << endl;
      if (prob <= ((1.+ConfLevel)/2.))
	break;
      --di;
    }
}
