#include "FitPainter.h"
#include "CommandParser.h"
#include "DrawLogo.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>

#include <stdlib.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TLine.h>
#include <math.h>
#include <TMath.h>

struct chi2type
{
  double chi2;
  double chi2_00;
  int dof;
};
typedef map<string, map<int, chi2type> > chi2listtype;

void FitPainter(vector<string> dirs)
{
  //make chi2 list
  chi2listtype chi2list;

  //read chi2
  for (vector<string>::iterator dit = dirs.begin(); dit != dirs.end(); dit++)
    {
      string fname = *dit + "/Results.txt";
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
	      chi2type ch;
	      ch.chi2 = chi2value;
	      ch.chi2_00 = -1;
	      ch.dof = dof;
	      chi2list["Total"][dit - dirs.begin()] = ch;
	    }

	  //read Correlated chi2
	  if (buffer == "Correlated")
	    {
	      iss >> dummy  >> chi2value;
	      chi2type ch;
	      ch.chi2 = chi2value;
	      ch.chi2_00 = -1;
	      ch.dof = 0;
	      chi2list["Correlated"][dit - dirs.begin()] = ch;
	    }

	  if (buffer != "Dataset")
	      continue;
	  iss >> index  >> chi2value  >> dof;

	  iss >> dummy;
	  name = dummy; //dataset name
	  while (iss >> dummy)
	    name += " " + dummy;
	  
	  chi2type ch;
	  ch.chi2 = chi2value;
	  ch.chi2_00 = -1;
	  ch.dof = dof;
	  chi2list[name][dit - dirs.begin()] = ch;
	}
      f.close();
    }
  //make chi2_00 list
  chi2listtype chi2list_00;

  //read chi2 without PDF uncertainties
  if (opts.chi2nopdf)
    for (vector<string>::iterator dit = dirs.begin(); dit != dirs.end(); dit++)
      {
	string fname_00 = *dit + "/Results_00.txt";
	ifstream f_00(fname_00.c_str());
	if (!f_00.good())
	  continue;

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
		chi2list["Total"][dit - dirs.begin()].chi2_00 = chi2value;
	      }

	    //read Correlated chi2
	    if (buffer == "Correlated")
	      {
		iss >> dummy  >> chi2value;
		chi2list["Correlated"][dit - dirs.begin()].chi2_00 = chi2value;
	      }

	    if (buffer != "Dataset")
	      continue;
	    iss >> index  >> chi2value  >> dof;

	    iss >> dummy;
	    name = dummy;
	    while (iss >> dummy)
	      name += " " + dummy;
	  
	    chi2list[name][dit - dirs.begin()].chi2_00 = chi2value;
	  }
	f_00.close();
      }


  int ndata = chi2list.size();

  float chi2[dirs.size()][ndata];
  float chi2_00[dirs.size()][ndata];
  int dof[dirs.size()][ndata];
  string dataname[ndata];
  string dirname[dirs.size()];
  for (int i = 0; i < dirs.size(); i++)
    for (int d = 0; d < ndata; d++)
      {
	chi2[i][d] = 0;
	chi2_00[i][d] = 0;
	dof[i][d] = 0;
      }

  for (int d = 0; d < dirs.size(); d++)
    {
      dirname[d] = opts.labels[d];
      while (dirname[d].find("/") != string::npos)
	dirname[d].replace(dirname[d].find("/"), 1, " ");
      while (dirname[d].find("_") != string::npos)
	dirname[d].replace(dirname[d].find("_"), 1, " ");
    }

  //loop on datasets
  int dat = 0;
  for (chi2listtype::iterator cit = chi2list.begin(); cit != chi2list.end(); cit++)
    {
      if ((*cit).first == "Correlated" || (*cit).first == "Total")
	continue;
      dataname[dat] = (*cit).first;
      //loop on directories
      int dir = 0;
      for (map<int, chi2type>::iterator dit = (*cit).second.begin(); dit != (*cit).second.end(); dit++)
	{
	  chi2[dir][dat] = dit->second.chi2;
	  chi2_00[dir][dat] = dit->second.chi2_00;
	  dof[dir][dat] = dit->second.dof;
	  dir++;
	}
      dat++;
    }

  dataname[dat] = "Correlated";
  //loop on directories
  int dir = 0;
  for (map<int, chi2type>::iterator dit = chi2list.find("Correlated")->second.begin(); dit != chi2list.find("Correlated")->second.end(); dit++)
    {
      chi2[dir][dat] = dit->second.chi2;
      chi2_00[dir][dat] = dit->second.chi2_00;
      dof[dir][dat] = dit->second.dof;
      dir++;
    }
  dat++;
  dataname[dat] = "Total";
  //loop on directories
  dir = 0;
  for (map<int, chi2type>::iterator dit = chi2list.find("Total")->second.begin(); dit != chi2list.find("Total")->second.end(); dit++)
    {
      chi2[dir][dat] = dit->second.chi2;
      chi2_00[dir][dat] = dit->second.chi2_00;
      dof[dir][dat] = dit->second.dof;
      dir++;
    }

  //write table
  FILE *ftab = fopen((opts.outdir + "chi2.tex").c_str(), "w");
  if (!ftab)
    {
      cout << "Cannot open chi2 tex file" << endl;
      return;
    }

  int points;
  if (dirs.size() >= 1)
    points = 14;
  if (dirs.size() >= 3)
    points = 12;
  if (dirs.size() >= 4)
    points = 11;
  if (dirs.size() >= 5)
    points = 9;

  float cm = 0.035277778 * points * 5.2;

  fprintf(ftab,"\\documentclass[%dpt]{report}\n", points);
  fprintf(ftab,"\\usepackage{extsizes}\n");
  fprintf(ftab,"\\usepackage[paperwidth=21cm,paperheight=21cm,left=0.2cm,right=0.2cm]{geometry}\n");
  fprintf(ftab,"\\usepackage[scaled]{helvet}\n");
  fprintf(ftab,"\\renewcommand\\familydefault{\\sfdefault}\n");
  fprintf(ftab,"\\pagestyle{empty}\n");

  fprintf(ftab,"\\usepackage{booktabs}\n");
  fprintf(ftab,"\\begin{document}\n");

  fprintf(ftab,"\\begin{table}\n");
  fprintf(ftab,"  \\begin{center}\n");
  fprintf(ftab,"    \\begin{tabular}{l");
  for (int i = 0; i < dirs.size(); i++)
    fprintf(ftab,"p{%.2fcm}", cm);
  fprintf(ftab,"}\n");
  fprintf(ftab,"      \\toprule\n");

  fprintf(ftab,"  Dataset  ");
  for (int i = 0; i < dirs.size(); i++)
    fprintf(ftab,"   & %s", dirname[i].c_str());
  fprintf(ftab,"  \\\\ \n");
  fprintf(ftab,"      \\midrule\n");

  for (int d = 0; d < ndata - 2; d++)
    {
      fprintf(ftab,"  %s ", dataname[d].c_str());
      for (int i = 0; i < dirs.size(); i++)
	if (dof[i][d] != 0 || chi2[i][d] != 0 )
	  {
	    if (chi2_00[i][d] >= 0)
	      fprintf(ftab,"& %s (%s) / %d", Round(chi2[i][d])[0].c_str(), Round(chi2_00[i][d])[0].c_str(), dof[i][d]);
	    else
	      fprintf(ftab,"& %s / %d", Round(chi2[i][d])[0].c_str(), dof[i][d]);
	  }
	else
	  fprintf(ftab,"& - ");
      fprintf(ftab,"  \\\\ \n");
    }

  fprintf(ftab,"  Correlated $\\chi^2$  ");
  for (int i = 0; i < dirs.size(); i++)
    if (chi2_00[i][ndata-2] >= 0)
      fprintf(ftab,"& %s (%s)", Round(chi2[i][ndata-2])[0].c_str(), Round(chi2_00[i][ndata-2])[0].c_str());
    else
      fprintf(ftab,"& %s", Round(chi2[i][ndata-2])[0].c_str());
  fprintf(ftab,"  \\\\ \n");

  bool pdfunc = false;
  fprintf(ftab,"      \\midrule\n");
  fprintf(ftab,"  Total $\\chi^2$ / dof  ");
  for (int i = 0; i < dirs.size(); i++)
    if (chi2_00[i][ndata-1] >= 0)
      {
	fprintf(ftab,"& %s (%s) / %d", Round(chi2[i][ndata-1])[0].c_str(), Round(chi2_00[i][ndata-1])[0].c_str(), dof[i][ndata-1]);
	pdfunc = true;
      }
    else
      fprintf(ftab,"& %s / %d", Round(chi2[i][ndata-1])[0].c_str(), dof[i][ndata-1]);
  fprintf(ftab,"  \\\\ \n");

  fprintf(ftab,"      \\midrule\n");
  fprintf(ftab,"  $\\chi^2$ p-value  ");
  for (int i = 0; i < dirs.size(); i++)
    fprintf(ftab,"& %.2f ", TMath::Prob(chi2[i][ndata-1], dof[i][ndata-1]));
  fprintf(ftab,"  \\\\ \n");


  fprintf(ftab,"      \\bottomrule\n");
  fprintf(ftab,"    \\end{tabular}\n");
  fprintf(ftab,"  \\end{center}\n");
  fprintf(ftab,"\\end{table}\n");
  if (pdfunc)
    fprintf(ftab,"$\\chi^2$ with (w/o) PDF uncertainties  \n");


  fprintf(ftab,"\n");
  fprintf(ftab,"\\end{document}\n");
  fclose(ftab);

  string command = "pdflatex -interaction=batchmode -output-directory=" + opts.outdir + " " + opts.outdir + "chi2.tex >/dev/null";
  //debug latex
  //string command = "pdflatex  -output-directory=" + opts.outdir + " " + opts.outdir + "chi2.tex";
  system(command.c_str());

  //cleanup
  string clean = "rm -f " 
    + opts.outdir + "chi2.aux " 
    + opts.outdir + "chi2.log " 
    + opts.outdir + "chi2.nav " 
    + opts.outdir + "chi2.out " 
    + opts.outdir + "chi2.snm " 
    + opts.outdir + "chi2.toc ";
  system(clean.c_str());

  return;
}
