#include "Par.h"

#include "Outdir.h"
#include "CommandParser.h"

#include <TMath.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

Par::Par(string dirname, string label)
{
  if (outdirs[label].IsMCreplica())
    {
      //loop on MC replica directories
      for (vector <string>::iterator it = outdirs[label].dirlist.begin(); it != outdirs[label].dirlist.end(); it++)
	{
	  string fname = (*it) + "/parsout_0";
	  ifstream f(fname.c_str());
	  if (!f.good()) //Something wrong, parameter file should be there
	    {
	      cout << "File " << fname << " not found" << endl;
	      exit(1);
	    }
	  string line;
	  int idx = 0;
	  getline(f, line);
	  while(!f.eof()) 
	    {
	      istringstream iss(line);
	      
	      int index;
	      string name;
	      double val, err;
	      iss >> index >> name >> val >> err;

	      if (val != 0 && err == 0)
		{
		  parlist[index].name = name;
		  parlist[index].value = val;
		  parlist[index].error_p = err;
		  parlist[index].error_m = err;
		}
	      else if (val != 0)
		{
		  parlist[index].name = name;
		  MCparams[index].push_back(val);
		}

	      getline(f, line);
	    }
	  f.close();
	}
      fitstatus = "MC-replica";
      if (outdirs[label].IsMedian())
	uncertainties = "median";
      else
	uncertainties = "mean";

      if (outdirs[label].Is68cl())
	uncertainties += "$\\pm$68cl";
      else if (outdirs[label].Is90cl())
	uncertainties += "$\\pm$90cl";
      else
	uncertainties += "$\\pm$rms";

      //loop on parameters
      for (map <int, vector <double> >::iterator it = MCparams.begin(); it != MCparams.end(); it++)
	{
	  int index = it->first;
	  vector <double> xi;
	  //loop on MC replica
	  for (vector <double>::iterator mcit = it->second.begin(); mcit != it->second.end(); mcit++)
	    xi.push_back(*mcit);

	  double val = 0;
	  double eplus = 0;
	  double eminus = 0;
	  if (outdirs[label].IsMedian())
	    val = Median(xi);
	  else
	    val = TMath::Mean(xi.begin(), xi.end());

	  if (outdirs[label].Is68cl())
	    {
	      if (outdirs[label].IsAsym())
		deltaasym(xi, val, eplus, eminus, cl(1));
	      else
		eplus = eminus = delta(xi, val, cl(1));
	    }
	  else if (outdirs[label].Is90cl())
	    {
	      if (outdirs[label].IsAsym())
		deltaasym(xi, val, eplus, eminus, 0.90);
	      else
		eplus = eminus = delta(xi, val, 0.90);
	    }
	  else
	    eminus = eplus = TMath::RMS(xi.begin(), xi.end());

	  parlist[index].value = val;
	  parlist[index].error_p = eplus;
	  parlist[index].error_m = eminus;
	}
    }
  else
    {
      string fname = dirname + "/parsout_0";
      ifstream f(fname.c_str());
      if (!f.good()) //may be this is not a fit
	{
	  cout << "File " << fname << " is empty (or io error)" << endl;
	  return;
	}
      
      string line;
      int idx = 0;
      getline(f, line);
      while(!f.eof()) 
	{
	  istringstream iss(line);
	  
	  int index;
	  string name;
	  double val, err;
	  iss >> index >> name >> val >> err;
	  
	  if (val != 0)
	    {
	      parlist[index].value = val;
	      parlist[index].error_p = err;
	      parlist[index].error_m = err;
	      parlist[index].name = name;
	    }
	  getline(f, line);
	}
      f.close();

      //Read fit status
      fitstatus = fitstat(dirname);
      uncertainties = "migrad-hesse";
    }
}

string fitstat(string dir)
{
  //Read fit status
  string status = "not-a-fit";
  string fname = dir + "/minuit.out.txt";
  ifstream ff(fname.c_str());
  if (!ff.good())
    return status;
  else
    status = "undefined";
  string line;
  while (getline(ff, line))
    {
      if (line.find("STATUS=OK") != string::npos)
	{
	  status = "converged";
	  break;
	}
      if (line.find("STATUS=FAILED") != string::npos)
	{
	  status = "failed";
	  break;
	}
      if (line.find("STATUS=NOT POSDEF") != string::npos)
	status = "pos-def-forced";
    }
  ff.close();
  return status;
}
