#include "Outdir.h"

#include "CommandParser.h"

#include <dirent.h>
#include <stdlib.h>
#include <algorithm>
#include <fstream>

map <string, Outdir> outdirs;

map <string, PdfData> pdfmap;
map <string, Data> datamap;
map <string, Chi2> chi2map;
map <string, Par> parmap;
map <string, Chi2scanData> chi2scanmap;
map <string, int> ndatamap;

//map <string, PdfData> pdfmapProfiled;

vector <TGraphAsymmErrors*> allgraphs;

//Constructor, load all the directory data
Outdir::Outdir(string dir) : dirname(dir), MCreplica(false), median(opts.median), cl68(opts.cl68), cl90(opts.cl90), asym(opts.asym), bands3(false), profiled(opts.profile), reweighted(opts.reweight), GK(opts.GKweight), BAY(opts.BAYweight), scale68(opts.scale68), RotatedSet(0), rotated(false), singlePDF(false),singlePDFlabel(0)
{
  //parse dirs for flags and labels
  string pattern = "";
  string origdirname = dirname;
  bool doprefix = true;
  while (dirname.find(":") != string::npos && doprefix)
    {
      doprefix = false;
      string prefix = dirname.substr(0, dirname.find(":"));
      if (prefix == "MC")
	{
	  doprefix = true;
	  pattern = dirname.substr(dirname.find(":")+1, dirname.size() - dirname.find(":")-1);
	  MCreplica = true;
	  dirname.erase(0, dirname.find(":")+1);
      }
      if (prefix == "68cl")
	{
	  doprefix = true;
	  if (cl90 == true)
	    {
	      cout << "Options 68cl and 90cl are mutually exclusive, cannot use both" << endl;
	      exit(1);
	    }
	  cl68 = true;
	  median = true;
	  dirname.erase(0, dirname.find(":")+1);
      }
      if (prefix ==  "90cl")
	{
	  doprefix = true;
	  if (cl68 == true)
	    {
	      cout << "Options 68cl and 90cl are mutually exclusive, cannot use both" << endl;
	      exit(1);
	    }
	  cl90 = true;
	  median = true;
	  dirname.erase(0, dirname.find(":")+1);
      }
      if (prefix == "median")
	{
	  doprefix = true;
	  median = true;
	  dirname.erase(0, dirname.find(":")+1);
	}
      if (prefix == "asym")
	{
	  doprefix = true;
	  asym = true;
	  dirname.erase(0, dirname.find(":")+1);
	}
      if (prefix == "profile")
	{
	  doprefix = true;
	  profiled = true;
	  dirname.erase(0, dirname.find(":")+1);
	}
      if (prefix == "reweight-GK")
	{
	  doprefix = true;
	  reweighted = true;
	  GK = true;
	  dirname.erase(0, dirname.find(":")+1);
	}
     if (prefix == "reweight-BAY")
	{
	  doprefix = true;
	  reweighted = true;
	  BAY = true;
	  dirname.erase(0, dirname.find(":")+1);
	}
      if (prefix == "rotate") 
	{
	  doprefix = true;
	  rotated = true;
	  dirname.erase(0, dirname.find(":")+1);
	  prefix = dirname.substr(0,dirname.find(":"));
	  RotatedSet = atoi(prefix.c_str());
	  // std::cout << RotatedSet << "\n";
	  dirname.erase(0, dirname.find(":")+1);
	}
      if (prefix == "set") 
	{
	  doprefix = true;
	  singlePDF = true;
	  dirname.erase(0, dirname.find(":")+1);
	  prefix = dirname.substr(0,dirname.find(":"));
	  singlePDFlabel = atoi(prefix.c_str());
	  // std::cout << RotatedSet << "\n";
	  dirname.erase(0, dirname.find(":")+1);
	}
      if (prefix == "scale68")
	{
	  doprefix = true;
	  scale68 = true;
	  dirname.erase(0, dirname.find(":")+1);
        }
      if (prefix == "3bands")
        {
          doprefix = true;
          bands3 = true;
          dirname.erase(0, dirname.find(":")+1);
        }
    }

  //now parse for the label
  if ((dirname.rfind(":") != string::npos))
    {
      label = dirname.substr(dirname.rfind(":")+1, dirname.size() - dirname.rfind(":") - 1);
      opts.labels.push_back(label);
      dirname.erase(dirname.rfind(":"), dirname.size());
      //if ( opts.profiled)
      //opts.labels.push_back(label+" profiled");	
    }
  else
    {
      label = origdirname;
      opts.labels.push_back(label);
      //if ( opts.profiled)
      //opts.labels.push_back(label+" profiled");	
    }

  //If plotting only one directory and no outdir name is provided, set the outdir to the current directory
  if (opts.dirs.size() == 1 && opts.outdir == "" && ! MCreplica)
    opts.outdir = dirname;

  if (MCreplica) //make array of subdirectories      
  {
      int NTotRep = 0;
      //if the given pattern is a directory, list all the subdirectories
      DIR *dir;
      dir = opendir (dirname.c_str());
      if (dir != NULL)
      {
          struct dirent *dp;
          while ((dp = readdir (dir)) != NULL)
          {
              string subdir = dp->d_name;
              DIR *sdir;
              sdir = opendir((dirname + "/" + subdir).c_str());
              if (sdir != NULL)
              {
                  // First check if directory contains something useful
                  ifstream infile((dirname + "/" + subdir + "/Results.txt").c_str());
                  if (infile.is_open()) {  
                      NTotRep += 1;
                      if ( opts.looseRepSelection ) {
                          // Good enough to try
                          dirlist.push_back(dirname + "/" + subdir + "/");
                      }
                  }
                  
                  if ( ! opts.looseRepSelection) {
                      // extensive fit status
                      string status = fitstat(dirname + "/" + subdir + "/") ;   
                      if (status == "converged" || status == "pos-def-forced" || opts.looseRepSelection ) {    
                          dirlist.push_back(dirname + "/" + subdir + "/");
                      }
                  }
                  
                  
                  closedir (sdir);
              }
          }
          closedir(dir);
      }
      else
      {
          //search one directory up, for all directories matching the given pattern:
          dirname = "./" + dirname;
          string pattern = dirname.substr(dirname.rfind("/") + 1,dirname.size() - dirname.rfind("/") -1);
          dirname.erase(dirname.rfind("/"),dirname.size() - dirname.rfind("/"));
          DIR *dir;
          dir = opendir((dirname).c_str());
          if (dir == NULL)
          {
              cout << "Error: invalid MC replic pattern, directory " << dirname << " does not exist" << endl;
              exit(1);
          }
          
          struct dirent *dp;
          while ((dp = readdir(dir)) != NULL)
          {
              string subdir = dp->d_name;
              if (subdir.find(pattern) != string::npos)
              {
                  DIR *sdir;
                  sdir = opendir((dirname + "/" + subdir).c_str());
                  if (sdir != NULL)
                  {
                      //check if Results.txt exists
                      ifstream infile((dirname + "/" + subdir + "/Results.txt").c_str());
                      if (infile.is_open())
                      {
                          NTotRep += 1;
                          if ( opts.looseRepSelection ) {
                              // Good enough to try
                              dirlist.push_back(dirname + "/" + subdir + "/");
                          }
                          else
                          {
                              // extensive fit status
                              string status = fitstat(dirname + "/" + subdir + "/") ;   
                              if (status == "converged" || status == "pos-def-forced" || opts.looseRepSelection ) {    
                                  dirlist.push_back(dirname + "/" + subdir + "/");
                              }
                          }
                          infile.close();
                      }
                      closedir (sdir);
                  }
              }
          }
          closedir(dir);
      }
      if (dirlist.size() == 0)
      {
          cout << "Error: No directories matchin pattern " << dirname << endl;
          exit(1);
      };
      cout << "Selected number of MC replica="  << dirlist.size() << " out of " << NTotRep << " directories with Results.txt file present  \n";
      
  } // end of MC replica
  
  else //This is a single directory, check if it exists
    {
      DIR *dir;
      dir = opendir (dirname.c_str());
      if (dir != NULL)
	closedir (dir);
      else
	{
	  cout << "Error: directory " << dirname << " not found" << endl;
	  exit(1);
	}
    }

  outdirs[label] = *this;

  //Load PDF data
  if (!opts.nopdfs)
    {
      PdfData pdf(dirname, label);
      pdfmap[label] = pdf;

      /*      if (opts.profiled) {
	PdfData Profiled(pdf,dirname,label);
	pdfmapProfiled[label] = Profiled;
	}*/
    }

  //Load datasets
  Data data(dirname, label);
  datamap[label] = data;

  //Load chi2
  Chi2 chi2(dirname, label);
  chi2map[label] = chi2;

  //Load parameters
  Par par(dirname, label);
  parmap[label] = par;

  //Load chi2scan
  Chi2scanData chi2scan(dirname);
  chi2scanmap[label] = chi2scan;
}

vector <float> q2list()
{
  vector <float> q2;
  for (map <string, PdfData>::iterator pit = pdfmap.begin(); pit != pdfmap.end(); pit++)
    for (map <float, Pdf>::iterator qit = pit->second.Central.begin(); qit != pit->second.Central.end(); qit++)
      q2.push_back(qit->first);

  sort (q2.begin(), q2.end());
  vector<float>::iterator it = unique (q2.begin(), q2.end());
  q2.resize(distance(q2.begin(), it));

  return q2;
}

vector <int> datalist()
{
  vector <int> dl;
  for (map <string, Data>::iterator dit = datamap.begin(); dit != datamap.end(); dit++) //loop on directories
    for (map <int, Dataset>::iterator iit = dit->second.datamap.begin(); iit != dit->second.datamap.end(); iit++) //loop on dataset index
      for (map <int, Subplot>::iterator sit = iit->second.subplots.begin(); sit != iit->second.subplots.end(); sit++) //loop on subplots
	dl.push_back(sit->first + iit->first * 100);

  sort (dl.begin(), dl.end());
  vector<int>::iterator it = unique (dl.begin(), dl.end());
  dl.resize(distance(dl.begin(), it));

  return dl;
}

vector <int> chi2list()
{
  vector <int> chi2;
  for (map <string, Chi2>::iterator cit = chi2map.begin(); cit != chi2map.end(); cit++) //loop on directories
    for (map <int, chi2type>::iterator it = cit->second.chi2list.begin(); it != cit->second.chi2list.end(); it++) //loop on dataset index
      chi2.push_back(it->first);

  sort (chi2.begin(), chi2.end());
  vector<int>::iterator it = unique (chi2.begin(), chi2.end());
  chi2.resize(distance(chi2.begin(), it));

  return chi2;
}

vector <int> parlist()
{
  vector <int> par;
  for (map <string, Par>::iterator dit = parmap.begin(); dit != parmap.end(); dit++) //loop on directories
    for (map <int, partype>::iterator it = dit->second.parlist.begin(); it != dit->second.parlist.end(); it++) //loop on parameters
      par.push_back(it->first);

  sort (par.begin(), par.end());
  vector<int>::iterator it = unique (par.begin(), par.end());
  par.resize(distance(par.begin(), it));

  return par;
}

int finddataindex(string name)
{
  int index = -1;

  for (map <string, Data>::iterator dit = datamap.begin(); dit != datamap.end() && index == -1; dit++) //loop on directories
    for (map <int, Dataset>::iterator iit = dit->second.datamap.begin(); iit != dit->second.datamap.end() && index == -1; iit++) //loop on dataset index
      if (iit->second.GetName() == name)
	index = iit->second.GetIndex();

  return index;
}
string finddataname(int index)
{
  string name = "";

  for (map <string, Data>::iterator dit = datamap.begin(); dit != datamap.end() && name == ""; dit++) //loop on directories
    for (map <int, Dataset>::iterator iit = dit->second.datamap.begin(); iit != dit->second.datamap.end() && name == ""; iit++) //loop on dataset index
      if (iit->second.GetIndex() == index)
	name = iit->second.GetName();

  return name;
}
string findparname(int index)
{
  string name = "";

  for (map <string, Par>::iterator dit = parmap.begin(); dit != parmap.end() && name == ""; dit++) //loop on directories
    for (map <int, partype>::iterator iit = dit->second.parlist.begin(); iit != dit->second.parlist.end() && name == ""; iit++) //loop on parameters index
      if (iit->first == index)
	name = iit->second.name;

  return name;
}
vector <string> chi2scanlist()
{
  vector <string> c2scanlist;
  for (map <string, Chi2scanData>::iterator dit = chi2scanmap.begin(); dit != chi2scanmap.end(); dit++) //loop on directories
    if (dit->second.chi2.size() > 0)
      c2scanlist.push_back(dit->second.label);

  sort (c2scanlist.begin(), c2scanlist.end());
  vector<string>::iterator it = unique (c2scanlist.begin(), c2scanlist.end());
  c2scanlist.resize(distance(c2scanlist.begin(), it));

  return c2scanlist;
}
