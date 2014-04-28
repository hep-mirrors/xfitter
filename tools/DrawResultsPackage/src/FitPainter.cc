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

struct chi2type
{
  double chi2;
  int dof;
  int dirid;
};
typedef map<string,vector<chi2type> > chi2listtype;

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
	      ch.dof = dof;
	      ch.dirid = dit - dirs.begin();
	      chi2list["Total"].push_back(ch);
	    }

	  //read Correlated chi2
	  if (buffer == "Correlated")
	    {
	      iss >> dummy  >> chi2value;
	      chi2type ch;
	      ch.chi2 = chi2value;
	      ch.dof = 0;
	      ch.dirid = dit - dirs.begin();
	      chi2list["Correlated"].push_back(ch);
	    }

	  if (buffer != "Dataset")
	      continue;
	  iss >> index  >> chi2value  >> dof;

	  iss >> dummy;
	  name = dummy;
	  while (iss >> dummy)
	    name += " " + dummy;
	  
	  chi2type ch;
	  ch.chi2 = chi2value;
	  ch.dof = dof;
	  ch.dirid = dit - dirs.begin();
	  chi2list[name].push_back(ch);
	}
      f.close();
    }
  int ndata = chi2list.size();

  float chi2[dirs.size()][ndata];
  int dof[dirs.size()][ndata];
  string dataname[ndata];
  string dirname[dirs.size()];
  for (int i = 0; i < dirs.size(); i++)
    for (int d = 0; d < ndata; d++)
      {
	chi2[i][d] = 0;
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
  int d = 0;
  for (chi2listtype::iterator cit = chi2list.begin(); cit != chi2list.end(); cit++)
    {
      if ((*cit).first == "Correlated" || (*cit).first == "Total")
	continue;
      dataname[d] = (*cit).first;
      //loop on directories
      for (vector<chi2type>::iterator dit = (*cit).second.begin(); dit != (*cit).second.end(); dit++)
	{
	  chi2[(*dit).dirid][d] = (*dit).chi2;
	  dof[(*dit).dirid][d] = (*dit).dof;
	}
      d++;
    }

  dataname[d] = "Correlated";
  for (vector<chi2type>::iterator dit = chi2list.find("Correlated")->second.begin(); dit != chi2list.find("Correlated")->second.end(); dit++)
    {
      chi2[(*dit).dirid][d] = (*dit).chi2;
      dof[(*dit).dirid][d] = (*dit).dof;
    }
  d++;
  dataname[d] = "Total";
  for (vector<chi2type>::iterator dit = chi2list.find("Total")->second.begin(); dit != chi2list.find("Total")->second.end(); dit++)
    {
      chi2[(*dit).dirid][d] = (*dit).chi2;
      dof[(*dit).dirid][d] = (*dit).dof;
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

  float cm = 0.035277778 * points * 5;

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

  for (int d = 0; d < ndata - 1; d++)
    {
      fprintf(ftab,"  %s ", dataname[d].c_str());
      for (int i = 0; i < dirs.size(); i++)
	if (dof[i][d] != 0 || chi2[i][d] != 0 )
	  fprintf(ftab,"& %.1f / %d", chi2[i][d], dof[i][d]);
	else
	  fprintf(ftab,"& - ");
      fprintf(ftab,"  \\\\ \n");
    }

  fprintf(ftab,"      \\midrule\n");
  fprintf(ftab,"  %s  ", dataname[ndata-1].c_str());
  for (int i = 0; i < dirs.size(); i++)
    fprintf(ftab,"& %.1f / %d", chi2[i][ndata-1], dof[i][ndata-1]);
  fprintf(ftab,"  \\\\ \n");
  fprintf(ftab,"      \\bottomrule\n");
  fprintf(ftab,"    \\end{tabular}\n");
  fprintf(ftab,"  \\end{center}\n");
  fprintf(ftab,"\\end{table}\n");

  fprintf(ftab,"\n");
  fprintf(ftab,"\\end{document}\n");
  fclose(ftab);

  string command = "pdflatex -interaction=batchmode -output-directory=" + opts.outdir + " " + opts.outdir + "chi2.tex >/dev/null";
  //debug latex
  //string command = "pdflatex  -output-directory=" + opts.outdir + " " + opts.outdir + "chi2.tex";
  system(command.c_str());

  return;
}
