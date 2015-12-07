#include "Chi2.h"

#include "Outdir.h"
#include "CommandParser.h"

#include <fstream>
#include <sstream>
#include <stdlib.h>

#include "FileOpener.h"

Chi2::Chi2(string dirname, string label)
{
  if (outdirs[label].IsMCreplica())
    {
      //here should make a cumulative plot of total chi2 of all replica
      //for (vector <string>::iterator it = outdirs[label].dirlist.begin(); it != outdirs[label].dirlist.end(); it++)
      //{}
      chi2tot.chi2 = 0;
      chi2tot.chi2_00 = 0;
      chi2tot.dof = 0;
      chi2corr.chi2 = 0;
      chi2corr.chi2_00 = 0;
      chi2corr.dof = 0;
      chi2log.chi2 = 0;
      chi2log.chi2_00 = 0;
      chi2log.dof = 0;
      return;
    }
 
  // string fname = dirname + "/Results.txt";
  // string fname0 = dirname + "/Results_0.txt";
  // ifstream f(fname.c_str());
  // if (!f) {
    // f.open(fname0.c_str());
  // }
  // if (!f.good())
    // {
      // cout << "File " << fname << " is empty (or io error)" << endl;
      // return;
    // }

  InFileOpener_t fo;
  fo.Add(dirname + "/Results.txt");
  fo.Add(dirname + "/Results_0.txt"); // --- Offset mode: central fit results
  // fo.Clear();
  if(fo.Open()) return;
  ifstream &f = fo.GetStream();
  string fname = fo.GetPath();

  //initialise chi2log to 0
  chi2log.chi2 = 0;
  chi2log.chi2_00 = 0;
  chi2log.dof = 0;

  //make chi2 list
  string line;
  string buffer;
  string dummy, name;
  int index, dof;
  float chi2value, chi2valuelog;
  string chi2read;
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
          chi2tot.chi2_00 = 0;
          chi2tot.dof = dof;
        }

      if (buffer == "Offset")
        {
          iss >> dummy;
          if(dummy == "corrected") {
            iss >> chi2value  >> dof;
            chi2tot.chi2 = chi2value;
            chi2tot.chi2_00 = 0;
            chi2tot.dof = dof;
          }
        }

      //read Correlated chi2
      if (buffer == "Correlated")
        {
          iss >> dummy  >> chi2value;
          chi2corr.chi2 = chi2value;
          chi2corr.chi2_00 = 0;
          chi2corr.dof = 0;
        }

      //read total log penalty chi2
      if (buffer == "Log")
        {
          iss >> dummy >> dummy  >> chi2value;
          chi2log.chi2 = chi2value;
          chi2log.chi2_00 = 0;
          chi2log.dof = 0;
        }

      if (buffer != "Dataset")
        continue;

      //read Partial chi2
      iss >> index  >> chi2read;
      if (chi2read.find("(") != string::npos) //the partial log term is given in brackets
        {
          if (chi2read.find("(") == (chi2read.size() - 1))
            {
              chi2read.erase(chi2read.find("("));
              chi2value = atof(chi2read.c_str());
              iss >> chi2read;
              if (chi2read.find(")") != string::npos)
                chi2read.erase(chi2read.find(")"));
              else
                {
                  cout << "Error parsing " << fname << endl;
                  return;
                }
              chi2valuelog = atof(chi2read.c_str());
              iss >> dof;
            }
          else
            {
              chi2valuelog = atof(chi2read.substr(chi2read.find("(") + 1, chi2read.size() - chi2read.find("(") - 1).c_str());
              chi2read.erase(chi2read.find("("), chi2read.size() - chi2read.find("("));
              chi2value = atof(chi2read.c_str());
              iss >> dof;
            }
        }
      else
        {
          chi2value = atof(chi2read.c_str());
          chi2valuelog = 0;
          iss >> dof;
        }

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
      if (dtindex == 0)
        {
          cout << "Error: could not find index for dataset " << name << endl;
          continue;
        }
      // chi2type ch;
      chi2list[dtindex].chi2 = chi2value;
      chi2list[dtindex].chi2_00 = 0;
      chi2list[dtindex].chi2_log = chi2valuelog;
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

          //read total log penalty chi2
          if (buffer == "Log")
            {
              iss >> dummy >> dummy  >> chi2value;
              chi2log.chi2_00 = chi2value;
            }

          if (buffer != "Dataset")
            continue;

          //read Partial chi2
          iss >> index  >> chi2read;
          if (chi2read.find("(") != string::npos)
            {
              chi2read.erase(chi2read.find("("));
              chi2value = atof(chi2read.c_str());
              iss >> chi2read;
              iss >> dof;
            }
          else
            {
              chi2value = atof(chi2read.c_str());
              iss >> dof;
            }

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
