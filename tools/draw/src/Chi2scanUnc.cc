#include "Chi2scanUnc.h"

#include "Outdir.h"
#include "CommandParser.h"
#include "DrawLogo.h"

/*
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
*/

bool Chi2scanUnc()
{
  vector <TCanvas*> cnvs;

  //This is a list of all chi2 scanned parameters, but the code currently work assuming there is only one parameter in the list of directories
  vector <string> list = chi2scanlist();
  if (list.size() == 0)
    return true;

  //List of all uncertainties
  int nunc = 6;

  string unc[6] = {"Statistical", "Systematic", "PDFs", "QCD scales", "Total (decomposed)", "Total (from fit)"};
  float unc_p[opts.labels.size()][nunc];
  float unc_m[opts.labels.size()][nunc];
  vector <string> dirname;
  dirname.resize(opts.labels.size());
  for (int i = 0; i < opts.labels.size(); i++)
    for (int d = 0; d < nunc; d++)
      {
	unc_p[i][d] = 0;
	unc_m[i][d] = 0;
      }

  for (int d = 0; d < opts.labels.size(); d++)
    {
      dirname[d] = opts.labels[d];
      while (dirname[d].find("/") != string::npos)
	dirname[d].replace(dirname[d].find("/"), 1, " ");
      while (dirname[d].find("_") != string::npos)
	dirname[d].replace(dirname[d].find("_"), 1, " ");
    }

  //loop on labels (directories)
  int dir = 0;
  for (vector<string>::iterator itl = opts.labels.begin(); itl != opts.labels.end(); itl++)
    {
      unc_p[dir][0] = chi2scanmap[*itl].statp;
      unc_m[dir][0] = chi2scanmap[*itl].statm;
      unc_p[dir][1] = chi2scanmap[*itl].systp;
      unc_m[dir][1] = chi2scanmap[*itl].systm;
      unc_p[dir][2] = chi2scanmap[*itl].pdfp;
      unc_m[dir][2] = chi2scanmap[*itl].pdfm;
      unc_p[dir][3] = chi2scanmap[*itl].scalep;
      unc_m[dir][3] = chi2scanmap[*itl].scalem;
      unc_p[dir][4] = chi2scanmap[*itl].totdecp;
      unc_m[dir][4] = chi2scanmap[*itl].totdecm;
      unc_p[dir][5] = chi2scanmap[*itl].totp;
      unc_m[dir][5] = chi2scanmap[*itl].totm;
      dir++;
    }


  //write table
  FILE *ftab = fopen((opts.outdir + "unc_summary.tex").c_str(), "w");
  if (!ftab)
    {
      cout << "Cannot open unc_summary tex file" << endl;
      return true;
    }

  int points;
  if (opts.labels.size() >= 1)
    points = 14;
  if (opts.labels.size() >= 3)
    points = 12;
  if (opts.labels.size() >= 4)
    points = 11;
  if (opts.labels.size() >= 5)
    points = 9;

  float cm = 0.035277778 * points * 9;

  float width = opts.labels.size() * cm + 9*points*0.035277778;
  float height = (nunc + 3)*points*0.035277778;

  fprintf(ftab,"\\documentclass[%dpt]{report}\n", points);
  fprintf(ftab,"\\usepackage{extsizes}\n");
  fprintf(ftab,"\\usepackage[paperwidth=21cm,paperheight=21cm,left=0.2cm,right=0.2cm,top=0.2cm,bottom=0.2cm]{geometry}\n");
  if (opts.font == "helvet")
    {
      fprintf(ftab,"\\usepackage[scaled]{helvet}\n");
      fprintf(ftab,"\\renewcommand\\familydefault{\\sfdefault}\n");
      fprintf(ftab,"\\usepackage{mathastext}\n");
    }
  else if (opts.font == "modernbright")
    fprintf(ftab,"\\usepackage{cmbright}\n");
  else if (opts.font == "palatino")
    fprintf(ftab,"\\usepackage{pxfonts}\n");
  fprintf(ftab,"\\pagestyle{empty}\n");

  fprintf(ftab,"\\usepackage{graphicx}\n");
  fprintf(ftab,"\\usepackage{booktabs}\n");
  fprintf(ftab,"\\usepackage[table]{xcolor}\n");
  fprintf(ftab,"\\definecolor{lightgray}{gray}{0.85}\n");
  fprintf(ftab,"\\usepackage{color}\n");
  fprintf(ftab,"\\begin{document}\n");

  fprintf(ftab,"\\begin{table}\n");
  fprintf(ftab,"  \\begin{center}\n");
  fprintf(ftab,"  \\rowcolors{2}{lightgray}{}\n");
  if (height >= width)
    fprintf(ftab,"  \\resizebox*{!}{\\textwidth}{\n");
  else
    fprintf(ftab,"  \\resizebox*{\\textwidth}{!}{\n");
  fprintf(ftab,"    \\begin{tabular}{l");
  for (int i = 0; i < opts.labels.size(); i++)
    fprintf(ftab,"p{%.2fcm}", cm);
  fprintf(ftab,"}\n");
  fprintf(ftab,"      \\toprule\n");
  //fprintf(ftab," Uncertainty  ");
  for (int i = 0; i < opts.labels.size(); i++)
    fprintf(ftab," & %s", dirname[i].c_str());
  fprintf(ftab,"  \\\\ \n");
  fprintf(ftab,"      \\midrule\n");

  fprintf(ftab,"  $%s$ ", chi2scanmap[*(opts.labels.begin())].label.c_str());
  for (vector<string>::iterator itl = opts.labels.begin(); itl != opts.labels.end(); itl++)
    //fprintf(ftab," &  %s ", Round(chi2scanmap[*itl].min, min(chi2scanmap[*itl].deltap,chi2scanmap[*itl].deltam)/10.)[0].c_str());
    fprintf(ftab," &  %s ", Round(chi2scanmap[*itl].min, min(chi2scanmap[*itl].deltap,chi2scanmap[*itl].deltam))[0].c_str());
  fprintf(ftab,"\\\\ \n");
  fprintf(ftab,"      \\midrule\n");

  for (int u = 0; u < nunc; u++)
    {
      if (u == nunc-1)
	fprintf(ftab,"      \\midrule\n");
      fprintf(ftab,"  %s ", unc[u].c_str());
      for (vector<string>::iterator itl = opts.labels.begin(); itl != opts.labels.end(); itl++)
	{
	  int l = itl-opts.labels.begin();
	  if (unc_p[l][u] != 0 && unc_m[l][u] != 0)
	    {
	      if (unc_p[l][u] == unc_m[l][u]) //symmetric case
		fprintf(ftab,"& $\\pm %s$", Round(unc_p[l][u], unc_p[l][nunc-1]/10.)[0].c_str());
	      else //asymmetric case
		fprintf(ftab,"& $+%s %s$", 
			//Round(unc_p[l][u], unc_p[l][nunc-1]/10.)[0].c_str(), 
			//Round(unc_m[l][u], unc_m[l][nunc-1]/10.)[0].c_str());
			Round(unc_p[l][u], unc_p[l][nunc-1])[0].c_str(), 
			Round(unc_m[l][u], unc_m[l][nunc-1])[0].c_str());
	    }
	  else
	    fprintf(ftab,"& - ");
	}
      fprintf(ftab,"  \\\\ \n");
    }

  fprintf(ftab,"      \\bottomrule\n");
  fprintf(ftab,"    \\end{tabular}\n");
  fprintf(ftab,"  }\n");
  fprintf(ftab,"  \\end{center}\n");
  fprintf(ftab,"\\end{table}\n");

  fprintf(ftab,"\n");
  fprintf(ftab,"\\end{document}\n");
  fclose(ftab);

  string command = "pdflatex -interaction=batchmode -output-directory=" + opts.outdir + " " + opts.outdir + "unc_summary.tex > /dev/null 2>&1 /dev/null"; //suppress output
  //debug latex
  //string command = "pdflatex  -output-directory=" + opts.outdir + " " + opts.outdir + "par.tex";

  bool latexcmd = system(command.c_str());
  if (latexcmd)
    {
      cout << "pdflatex error in making: unc_summary.pdf" << endl;
      cout << "check the error with:" << endl;
      cout << command << endl;
    }

  //cleanup
  string clean = "rm -f " 
    + opts.outdir + "unc_summary.aux " 
    + opts.outdir + "unc_summary.log " 
    + opts.outdir + "unc_summary.nav " 
    + opts.outdir + "unc_summary.out " 
    + opts.outdir + "unc_summary.snm " 
    + opts.outdir + "unc_summary.toc ";
  if (system(clean.c_str())) {
    cout << "Error removing temporary files\n";
  }

  //return latexcmd;
  return false;
}
