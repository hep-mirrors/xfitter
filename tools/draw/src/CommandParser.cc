#include "CommandParser.h"

#include <TStyle.h>

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
  scale68(false),
  profile(false),
  reweight(false),
  GKweight(false),
  BAYweight(false),
  bw(false),
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
  plotsperpage(2),
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
  noupband(false),
  errbandcol(kYellow),
  points(false),
  theorylabel("Theory"),
  onlytheory(false),
  threlerr(false),
  ratiototheory(false),
  diff(false),
  twopanels(false),
  threepanels(false),
  multitheory(false),
  nothshifts(false),
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
  logpenalty(false),
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

  lstyl[0] = 1;
  lstyl[1] = 2;
  lstyl[2] = 3;
  lstyl[3] = 4;
  lstyl[4] = 5;
  lstyl[5] = 6;

  // tight MC replica selection by default:
  looseRepSelection = false;
  
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
	else if (*it == "--largetext")
	  {
	    txtsize = 0.05;
	    lmarg = 0.18;
	    bmarg = 0.13;
	    offset = 1.6;
	  }
	else if (*it == "--bw")
	  bw = true;
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
	else if (*it == "--partial-log-penalty")
	  logpenalty = true;
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
	  {
	    cms = true;
	    drawlogo = false;
	  }
	else if (*it == "--cms-preliminary")
	  {
	    cmspreliminary = true;
	    drawlogo = false;
	  }
	else if (*it == "--atlas")
	  {
	    atlas = true;
	    drawlogo = false;
	  }
	else if (*it == "--atlas-internal")
	  {
	    atlasinternal = true;
	    drawlogo = false;
	  }
	else if (*it == "--atlas-preliminary")
	  {
	    atlaspreliminary = true;
	    drawlogo = false;
	  }
	else if (*it == "--cdfii-preliminary")
	  {
	    cdfiipreliminary = true;
	    drawlogo = false;
	  }
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
	else if (*it == "--scale68")
	  scale68 = true;
	else if (*it == "--profile") {
	  dobands = true;
	  profile = true;
	}
	else if (*it == "--reweight-BAY") {
	  dobands = true;
	  reweight = true;
	  BAYweight = true;
	}
	else if (*it == "--reweight-GK") {
	  dobands = true;
	  reweight = true;
	  GKweight = true;
	}
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
	else if (*it == "--plots-per-page")
	  {
	    plotsperpage = atoi((*(it+1)).c_str());
	    allargs.erase(it+1);
	  }
	else if (*it == "--loose-mc-replica-selection")
	  looseRepSelection = true;
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
		col[1] = kOrange;
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
	    else if (pattern == 4)
	      {
		col[0] = kBlue + 1;
		col[1] = kAzure - 9;
		col[2] = kAzure + 3;
		col[3] = kAzure + 4;
		col[4] = kAzure + 5;
		col[5] = kAzure + 6;
	      }
	    else if (pattern == 5)
	      {
		col[0] = kOrange + 7;
		col[1] = kYellow;
		col[2] = kRed - 1;
		col[3] = kOrange + 4;
		col[4] = kOrange + 2;
		col[5] = kOrange -1;
	      }
	    else if (pattern == 6)
	      {
		col[0] = kGreen + 2;
		col[1] = kSpring - 9;
		col[2] = kGreen + 1;
		col[3] = kSpring + 4;
		col[4] = kSpring + 2;
		col[5] = kSpring + 7;
	      }
	    else if (pattern == 7)
	      {
		col[0] = kRed - 2;
		col[1] = kAzure - 9;
		col[2] = kSpring - 9;
		col[3] = kOrange + 7;
		col[4] = kSpring + 2;
		col[5] = kSpring + 7;
	      }
	    else if (pattern == 8)
	      {
		col[0] = kRed - 2;
		col[1] = kBlue + 1;
		col[2] = kGreen + 1;
		col[3] = kOrange + 7;
		col[4] = kSpring + 2;
		col[5] = kSpring + 7;
	      }


	      allargs.erase(it+1);
	  }
	else if (*it == "--therr")
	  therr = true;
	else if (*it == "--noupband")
	  noupband = true;
	else if (*it == "--greenband")
	  errbandcol = kGreen - 3;
	else if (*it == "--blueband")
	  errbandcol = kAzure - 9;
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
	else if (*it == "--theory-rel-errors")
	  {
	    onlytheory = true;
	    ratiototheory = true;
	    threlerr = true;
	  }
	else if (*it == "--ratio-to-theory")
	  ratiototheory = true;
	else if (*it == "--diff")
	  diff = true;
	else if (*it == "--2panels")
	  twopanels = true;
	else if (*it == "--3panels")
	  threepanels = true;
	else if (*it == "--multitheory")
	  multitheory = true;
	else if (*it == "--nothshifts")
	  nothshifts = true;
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
vector<string> Round(double value, double error, bool sign)
{
  vector <string> result;

  int decimal = 0;

  //If no error, value is rounded to two significant digits
  if (error == 0)
    error = value;

  if (error != 0)
    decimal = -log10(fabs(error)) + 2;
  decimal = max(0, decimal);

  char Dec[10];
  sprintf (Dec, "%d", decimal);
  string D = Dec;

  char Numb[50];
  if (sign)
    sprintf (Numb, ((string)"%+." + D + "f").c_str(), value);
  else
    sprintf (Numb, ((string)"%." + D + "f").c_str(), value);
  result.push_back(Numb);
  sprintf (Numb, ((string)"%." + D + "f").c_str(), error);
  result.push_back(Numb);
  
  return result;
}
