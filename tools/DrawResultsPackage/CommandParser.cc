#include "CommandParser.h"
#include <stdlib.h>
#include <iostream>
#include <TH1F.h>

float txtsize = 0.043;
float offset = 1.5;
float lmarg = 0.15;
//float bmarg = 0.2;

CommandParser::CommandParser(int argc, char **argv):
  dobands(false),
  asymbands(false),
  logx(true),
  filledbands(false),
  rmin(0),
  rmax(2),
  xmin(-1),
  xmax(-1),
  pdf(false),
  splitplots(false),
  ext("eps"),
  resolution(1200),
  therr(false),
  points(false),
  theorylabel("Theory"),
  threepanels(false),
  version(true),
  drawlogo(true),
  nodata(false),
  nopdfs(false),
  cms(false),
  cmspreliminary(false),
  atlas(false),
  atlaspreliminary(false),
  atlasinternal(false),
  cdfiipreliminary(false),
  outdir("")
{

  //initialise colors and styles
  colors[0] = kRed;
  colors[1] = kBlue + 2;
  colors[2] = kGreen + 3;
  colors[3] = kMagenta;
  colors[4] = kAzure + 1;
  colors[5] = kOrange + 7;

  styles[0] = 3004;
  styles[1] = 3005;
  styles[2] = 3006;
  styles[3] = 3007;
  styles[4] = 3016;
  styles[5] = 3020;

  markers[0] = 24;
  markers[1] = 25;
  markers[2] = 26;
  markers[3] = 32;
  markers[4] = 31;
  markers[5] = 27;

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
	else if (*it == "--pdf")
	  pdf = true;
	else if (*it == "--lowres")
	  resolution = 600;
	else if (*it == "--highres")
	  resolution = 2400;
	else if (*it == "--no-version")
	  version = false;
	else if (*it == "--no-logo")
	  drawlogo = false;
	else if (*it == "--no-data")
	  nodata = true;
	else if (*it == "--no-pdfs")
	  nopdfs = true;
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
	else if (*it == "--asymbands")
	  {
	    dobands = true;
	    asymbands = true;
	  }
	else if (*it == "--no-logx")
	  logx = false;
	else if (*it == "--outdir")
	  {
	    outdir = *(it+1);
	    allargs.erase(it+1);
	  }
	else if (*it == "--splitplots")
	  splitplots = true;
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
		colors[0] = kBlue;
		colors[1] = kYellow;
		colors[2] = kGreen;
		colors[3] = kRed;
		colors[4] = kMagenta;
		colors[5] = kCyan;
	      }
	    else if (pattern == 2)
	      {
		colors[0] = kBlue;
		colors[1] = kRed;
		colors[2] = kYellow;
		colors[3] = kOrange;
		colors[4] = kMagenta;
		colors[5] = kCyan;
	      }
	    else if (pattern == 3)
	      {
		colors[0] = kBlue;
		colors[1] = kMagenta;
		colors[2] = kCyan;
		colors[3] = kRed;
		colors[4] = kGreen + 2;
		colors[5] = kYellow;
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

  //parse dirs for labels
  for (vector<string>::iterator it = dirs.begin(); it != dirs.end(); it++)
    if ((*it).find(":") != string::npos)
      {
	labels.push_back((*it).substr((*it).find(":")+1, (*it).size() - (*it).find(":") - 1));
	(*it).erase((*it).find(":"), (*it).size());
      }
    else
      labels.push_back((*it));

  if (outdir == "")
    if (dirs.size() == 1) {outdir = dirs[0];}
    else {outdir = "plots/";}

  if (outdir.rfind("/") != outdir.size() - 1)
    outdir.append("/");
  
  cout << endl;
  cout << "plots are stored in: " << outdir << endl;
}

CommandParser opts;
