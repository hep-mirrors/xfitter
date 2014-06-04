#include "Chi2.h"

#include "Outdir.h"
#include "CommandParser.h"

#include <fstream>
#include <sstream>

Chi2::Chi2(string dirname, string label)
{
for (vector <string>::iterator it = outdirs[label].dirlist.begin(); it != outdirs[label].dirlist.end(); it++)

  if (outdirs[label].IsMCreplica())
    {
      //here should make a cumulative plot of total chi2 of all replica
      chi2tot.chi2 = 0;
      chi2tot.chi2_00 = -1;
      chi2tot.dof = 0;
      chi2corr.chi2 = 0;
      chi2corr.chi2_00 = -1;
      chi2corr.dof = 0;
      return;
    }

  string fname = dirname + "/Results.txt";
  ifstream f(fname.c_str());
  if (!f.good())
    {
      cout << "File " << fname << " is empty (or io error)" << endl;
      return;
    }

  //make chi2 list
  string line;
  string buffer;
  string dummy, name;
  int index, dof;
  float chi2value;
  while (getline(f, line))
    {
      buffer = "";
      istringstream iss(line);
      iss >> buffer; 
      
      //read total chi2
      if (buffer == "After")
	{
	  iss >> dummy  >> chi2value  >> dof;
	  chi2tot.chi2 = chi2value;
	  chi2tot.chi2_00 = -1;
	  chi2tot.dof = dof;
	}

      //read Correlated chi2
      if (buffer == "Correlated")
	{
	  iss >> dummy  >> chi2value;
	  chi2corr.chi2 = chi2value;
	  chi2corr.chi2_00 = -1;
	  chi2corr.dof = 0;
	}

      if (buffer != "Dataset")
	continue;
      iss >> index  >> chi2value  >> dof;

      //read Partial chi2
      iss >> dummy;
      if (line.find(dummy) == string::npos)
	{
	  cout << "Error parsing " << fname << endl;
	  return;
	}
      name = line.substr(line.find(dummy), line.size() - line.find(dummy));
      string trimmed = name;
      trimmed.erase(trimmed.find_last_not_of(" ")+1, string::npos);
      name = trimmed;
      int dtindex = finddataindex(name);
      if (dtindex == -1)
	{
	  cout << "Error: could not find index for dataset " << name << endl;
	  continue;
	}
      chi2type ch;
      chi2list[dtindex].chi2 = chi2value;
      chi2list[dtindex].chi2_00 = -1;
      chi2list[dtindex].dof = dof;
    }
  f.close();

  //read chi2 without PDF uncertainties
  string fname_00 = dirname + "/Results_00.txt";
  ifstream f_00(fname_00.c_str());
  if (opts.chi2nopdf && f_00.good())
    {
      //make chi2_00 list
      string line;
      string buffer;
      string dummy, name;
      int index, dof;
      float chi2value;
      while (getline(f_00, line))
	{
	  buffer = "";
	  istringstream iss(line);
	  iss >> buffer; 

	  //read total chi2
	  if (buffer == "After")
	    {
	      iss >> dummy  >> chi2value  >> dof;
	      chi2tot.chi2_00 = chi2value;
	    }

	  //read Correlated chi2
	  if (buffer == "Correlated")
	    {
	      iss >> dummy  >> chi2value;
	      chi2corr.chi2_00 = chi2value;
	    }

	  if (buffer != "Dataset")
	    continue;
	  iss >> index  >> chi2value  >> dof;

	  //read Partial chi2
	  iss >> dummy;
	  name = line.substr(line.find(dummy), line.size() - line.find(dummy));
	  string trimmed = name;
	  trimmed.erase(trimmed.find_last_not_of(" ")+1, string::npos);
	  name = trimmed;
	  int dtindex = finddataindex(name);
	  if (dtindex == -1)
	    {
	      cout << "Error: could not find index for dataset " << name << endl;
	      continue;
	    }
	  chi2list[dtindex].chi2_00 = chi2value;
	}
      f_00.close();
    }
}
