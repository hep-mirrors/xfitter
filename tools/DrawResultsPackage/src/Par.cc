#include "Par.h"

#include <iostream>
#include <fstream>
#include <sstream>

Par::Par(string dirname, string label)
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
	  parlist[index].error = err;
	  parlist[index].name = name;
	}
      getline(f, line);
    }
  f.close();

  //Read fit status
  fitstatus = "not-a-fit";
  fname = dirname + "/minuit.out.txt";
  ifstream ff(fname.c_str());
  if (!ff.good())
    return;
  else
    fitstatus = "converged";
  while (getline(ff, line))
    if (line.find("MATRIX FORCED POS-DEF") != string::npos)
      fitstatus = "pos-def-forced";
  ff.close();
}
