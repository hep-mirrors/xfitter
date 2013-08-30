#include "CommandParser.h"
#include <stdlib.h>
#include <iostream>
#include <TH1F.h>


CommandParser::CommandParser(int argc, char **argv):
  dobands(false),
  splitplots(false),
  filledbands(false),
  outdir(""),
  rmin(0),
  rmax(2)
{

  //initialise colors and styles
  colors[0] = kRed;
  colors[1] = kBlue;
  colors[2] = kMagenta;
  colors[3] = kGreen + 2;
  colors[4] = kCyan;
  colors[5] = kYellow;

  styles[0] = 3004;
  styles[1] = 3005;
  styles[2] = 3006;
  styles[3] = 3007;
  styles[4] = 3016;
  styles[5] = 3020;

  //read all command line arguments
  vector <string> allargs;
  for (int iar = 1; iar < argc; iar++)
    allargs.push_back(argv[iar]);

  //search for options
  for (vector<string>::iterator it = allargs.begin(); it != allargs.end(); it++)
    if ((*it).find("--") == 0)
      {
	if (*it == "--help")
	  {
	    cout << endl;
	    cout << "program usage:" << endl;
	    cout << argv[0] << " [options] dir1[:label1] [dir2:[label2]] [...]" << endl;
	    cout << endl;
	    cout << "options:" << endl;
	    cout << "\t --bands" << endl;
	    cout << "\t --outdir <output directory>" << endl;
	    cout << "\t --splitplots" << endl;
	    cout << "\t --filledbands" << endl;
	    cout << "\t --ratiorange <min:max>" << endl;
	    cout << "\t --colorpattern <1-3>" << endl;
	    cout << endl;
	    cout << " dir1 is used as reference for ratio plots" << endl;
	    exit(0);
	  }
	else if (*it == "--bands")
	  dobands = true;
	else if (*it == "--outdir")
	  {
	    outdir = *(it+1);
	    allargs.erase(it+1);
	  }
	else if (*it == "--splitplots")
	  splitplots = true;
	else if (*it == "--filledbands")
	  filledbands = true;
	else if (*it == "--ratiorange")
	  {
	    rmin = atof((*(it+1)).substr(0, (*(it+1)).find(":")).c_str());
	    cout << (*(it+1)).substr(0, (*(it+1)).find(":")) << endl;
	    rmax = atof((*(it+1)).substr((*(it+1)).find(":") + 1, (*(it+1)).size() - (*(it+1)).find(":") - 1).c_str());
	    cout << (*(it+1)).substr((*(it+1)).find(":") + 1, (*(it+1)).size() - (*(it+1)).find(":") - 1) << endl;
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
	else
	  {
	    cout << "Invalid option " << *it << endl;
	    exit(0);
	  }
	allargs.erase(it);
	it = allargs.begin();
      }
  
  for (vector<string>::iterator it = allargs.begin(); it != allargs.end(); it++)
    dirs.push_back(*it);

  if (dirs.size() == 0)
    {
      //------------------------
      cout << endl;
      cout << "program usage:" << endl;
      cout << argv[0] << " [options] dir1[:label1] [dir2:[label2]] [...]" << endl;
      cout << endl;
      cout << " Specify at least one directory" << endl;
      //------------------------
      exit(0);
    }



  if (dirs.size() > 6)
    {
      //------------------------
      cout << "program usage:" << endl;
      cout << argv[0] << " [options] dir1[:label1] [dir2:[label2]] [...]" << endl;
      cout << endl;
      cout << " maximum number of directories is 6" << endl;
      //------------------------
      exit(0);
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
