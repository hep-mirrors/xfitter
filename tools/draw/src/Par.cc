#include "Par.h"

#include "Outdir.h"
#include "CommandParser.h"
#include "pdferrors.h"

#include <TMath.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include "FileOpener.h"
#include "FTNFitPars.h"

std::vector<std::pair<double,double>> pumplin_errors(const std::string& dir);

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
            val = median(xi);
          else
            val = mean(xi);

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
            eminus = eplus = rms(xi);

          parlist[index].value = val;
          parlist[index].error_p = eplus;
          parlist[index].error_m = eminus;
        }
    }
  else
    {
      // string fname = dirname + "/parsout_0";
      // ifstream f(fname.c_str());
      // if (!f.good()) //may be this is not a fit
        // {
          // cout << "File " << fname << " is empty (or io error) (in par.cc)" << endl;
          // return;
        // }
        
      InFileOpener_t fo;
      fo.Add(dirname + "/parsout_0");
      fo.Add(dirname + "/MI_saved_final.txt");
      // fo.Clear();
      if(fo.Open()) return;
      
      if(fo.GetIndex() == 1) {
        // --- 'MI_saved_final.txt' opened
        // --- which means that Offset_Finalize has been called 
        fitstatus = "converged"; // ws: who needs this?
        // -- read parlist
        fo.Close();
        FTNFitPars_t mfp;
        mfp.Read(fo.GetPath().c_str(), "p");  
        vector<int> vi = mfp.GetVarIndices();
        for(vector<int>::iterator it = vi.begin(); it != vi.end(); it++) {
          // int index = *it;
          int index = mfp.UID(*it);
          parlist[index].value = mfp.Value(*it);
          parlist[index].error_p =
          parlist[index].error_m = mfp.Error(*it);
          parlist[index].name = mfp.Name(*it);
          // cout << index <<": "<< parlist[index].value << endl;
        }
      }
      else {
	// read pumplin errors:
	std::vector<std::pair<double,double>> eJon = pumplin_errors(dirname);

        ifstream &f = fo.GetStream();
        string line;
        int idx = 0;
        getline(f, line);
	int internalIndex = 0;
        while(!f.eof()) 
          {
            istringstream iss(line);
            
            int index; //--- Minuit parameter index
            string name;
            double val, err;
            iss >> index >> name >> val >> err;
            
            if (val != 0)
              {
                parlist[index].value = val;
                parlist[index].error_p = err;
                parlist[index].error_m = err;
                parlist[index].name = name;
		if ( (err !=0) and (eJon.size()>0) and (internalIndex<eJon.size()) ) {
		  // update error
		  parlist[index].error_p = eJon[internalIndex].first;
		  parlist[index].error_m = eJon[internalIndex].second;
		  internalIndex++;
		}
              }
            getline(f, line);
          }
        f.close();
        //Read fit status
        fitstatus = fitstat(dirname);
      }
      uncertainties = hessestat(dirname);
    }
}

string fitstat(string dir)
{
  //Read fit status
  string status = "not-a-fit";

  // First check the Status.out file:
  ifstream fOk((dir+"/Status.out").c_str());
  if (!fOk.good()) { 
    return status;
  }
  else {
    string line;
    while (getline(fOk,line)) {
      if (line.find("OK") !=  string::npos) 	{
	status = "fit OK";
      }
      else {
	status = "not-a-fit";
	return status;
      }
    }
  }

  string fname = dir + "/minuit.out.txt";
  ifstream ff(fname.c_str());
  if (!ff.good())
    return status;
  else
    status = "undefined";
  string line;
  while (getline(ff, line))
    {
      if (line.find("STATUS=CONVERGED") != string::npos)
        {
          status = "converged";
          break;
        }
      if (line.find("STATUS=FAILED") != string::npos)
          status = "failed";
    }
  ff.close();
  return status;
}

string hessestat(string dir)
{
  //Read fit status
  string status = "not-a-fit";

  // First check the Status.out file:
  ifstream fOk((dir+"/Status.out").c_str());
  if (!fOk.good()) { 
    return status;
  }
  else {
    string line;
    while (getline(fOk,line)) {
      if (line.find("OK") !=  string::npos) 	{
	status = "fit OK";
      }
      else {
	status = "not-a-fit";
	return status;
      }
    }
  }

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
	status = "migrad-hesse";  //--- a text printed by ParPainter
      if (line.find("STATUS=NOT POSDEF") != string::npos)
        status = "pos-def-forced";
      if (line.find("ERROR MATRIX FROM ITERATE") != string::npos)
	status = "iterate";
    }
  ff.close();
  return status;
}

std::vector<std::pair<double,double>> pumplin_errors(const std::string& dir) {
  std::vector<std::pair<double,double>> pars;

  string fname = dir + "/minuit.out.txt";
  ifstream ff(fname.c_str());
  if (!ff.good())
    return pars;
  // find "MYSTUFF" statement
  string line;
  bool gotMyStuff;
  while (getline(ff, line)) {
    if (line.find("mystuff called with argument  1000") != string::npos) {
      gotMyStuff = true;
      break;
    }
  };
  if (! gotMyStuff)
    return pars;
  while (getline(ff, line)) {
    double errSym;
    if (line.find("mystuff called with argument  2000") != string::npos) 
      return pars;  // stop
    if ( (line.find("for DeltaChisq") != string::npos) and (line.find(" at ") != string::npos)) {
      int pos1 = line.find("+/-")+3;
      int pos2 = line.find("for DeltaChisq")-1;
      istringstream iss(line.substr(pos1,pos2));
      iss >> errSym;
    }
    if ((line.find("Steps in ") != string::npos)) {
      int pos1 = line.find("=")+1;
      int pos2 = line.find(",");
      int pos3 = line.find("(")-1;
      double deltaUp,deltaDown;
      istringstream iss1(line.substr(pos1,pos2-1));
      iss1 >> deltaUp;
      istringstream iss2(line.substr(pos2+1,pos3));
      iss2 >> deltaDown;
      pars.push_back(std::pair<double,double>(sqrt(deltaUp)*errSym,sqrt(deltaDown)*errSym));
    }
  }
  return pars;
}
