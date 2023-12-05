#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TError.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>

#include "Outdir.h"
#include "CommandParser.h"
#include "PdfsPainter.h"
#include "DataPainter.h"
#include "ShiftPainter.h"
#include "FitPainter.h"
#include "ParPainter.h"
#include "Chi2scanPainter.h"
#include "Chi2scanUnc.h"
#include "Chi2scanGauss.h"

using namespace std;

int main(int argc, char **argv)
{
  //--------------------------------------------------
  //parse command line arguments
  opts = CommandParser(argc, argv);

  gErrorIgnoreLevel=1001;

  //read output directories
  for (vector<string>::iterator itd = opts.dirs.begin(); itd != opts.dirs.end(); itd++)
    Outdir out((*itd).c_str());

  //check there are no repetion in labels to avoid same name in root TH1
  for (vector<string>::iterator it1 = opts.labels.begin(); it1 != opts.labels.end(); it1++)
    for (vector<string>::iterator it2 = it1+1; it2 != opts.labels.end(); it2++)
      if (*it1 == *it2)
        {
          cout << endl;
          cout << "Error: label (or directory) " << *it1 << " can appear only once in labels list" << endl;
          cout << "Specify different labels" << endl;
          cout << endl;
          exit(-1);
        }

  //Associate colors and styles to labels
  for (vector<string>::iterator itl = opts.labels.begin(); itl != opts.labels.end(); itl++)
    {
      opts.colors[*itl]  = opts.col[itl-opts.labels.begin()];
      opts.styles[*itl]  = opts.styl[itl-opts.labels.begin()];
      opts.lstyles[*itl]  = opts.lstyl[itl-opts.labels.begin()];
      opts.markers[*itl] = opts.mark[itl-opts.labels.begin()];
    }

  //Set default out directory
  if (opts.outdir == "")
    opts.outdir = "plots/";

  if (opts.outdir.rfind("/") != opts.outdir.size() - 1)
    opts.outdir.append("/");

  //--------------------------------------------------
  //Pdf plots
  vector <TCanvas*> pdfscanvaslist, pdfscanvasratiolist;
  if (! opts.nopdfs)
    {
      //loop on Q2 bins
      vector <float> ql = q2list();
      for (vector<float>::iterator qit = ql.begin(); qit != ql.end(); qit++)
        {
          //loop on pdf types
          for (vector <pdftype>::iterator pit = pdfs.begin(); pit != pdfs.end(); pit++)
            {
              vector <TCanvas*> pdfcnv = PdfsPainter(*qit, *pit);
              pdfscanvaslist.push_back(pdfcnv[0]);
              if (pdfcnv.size() > 1)
                pdfscanvasratiolist.push_back(pdfcnv[1]);
            }

          if (!opts.q2all)
            break;
        }
    }

  //--------------------------------------------------
  //Data plots
  gStyle->SetEndErrorSize(4);
  vector <TCanvas*> datapullscanvaslist;
  if (! opts.nodata)
    {
      //loop on datasets
      vector <int> dl = datalist();
      for (vector<int>::iterator dit = dl.begin(); dit != dl.end(); dit++)
        {
          //extract dataset index and subplot index
          int dataindex = (int)(*dit) / 100;
          int subplotindex = *dit - dataindex * 100;
          TCanvas *dataplot = DataPainter(dataindex, subplotindex);
          if (dataplot != 0)
            datapullscanvaslist.push_back(dataplot);
        }
    }

  //--------------------------------------------------
  //Shift plots
  vector <TCanvas*> shiftcanvaslist;
  if (! opts.noshifts)
    shiftcanvaslist = ShiftPainter(opts.dirs);

  //--------------------------------------------------
  //Chi2scan plots
  vector <TCanvas*> chi2scancanvaslist;
  //  if (! opts.nochi2scan)
  chi2scancanvaslist = Chi2scanPainter();

  //--------------------------------------------------
  //Chi2gauss plots
  vector <TCanvas*> chi2scangausscanvaslist;
  //  if (! opts.nochi2scan)
  chi2scangausscanvaslist = Chi2scanGauss();
  
  //--------------------------------------------------
  //Create output directory
  system(((string)"mkdir -p " + opts.outdir).c_str());

  //--------------------------------------------------
  //Fit and parameters results plots
  bool chi2tab = true;
  bool partab = true;
  bool uncsummary = true;
  if (!opts.notables)
    {
      chi2tab = FitPainter();
      partab = ParPainter();
      uncsummary = Chi2scanUnc();
    }

  //--------------------------------------------------
  //Save plots
  int pgn = 0;
  char pgnum[15];
  int plotsperpage = 2;
  gStyle->SetPaperSize(opts.pagewidth*float(opts.plotsperpage)/2., opts.pagewidth*float(opts.plotsperpage)/2.);
  vector <TCanvas*>::iterator it = pdfscanvaslist.begin();
  for (vector <TCanvas*>::iterator it = pdfscanvaslist.begin(); it != pdfscanvaslist.end();)
    {
      char numb[25];
      sprintf(numb, "%ld", it - pdfscanvaslist.begin());
      TCanvas * pagecnv = new TCanvas(numb, "", 0, 0, opts.resolution * opts.plotsperpage, opts.resolution * opts.plotsperpage);
      pagecnv->Divide(opts.plotsperpage, opts.plotsperpage);
      for (int i = 1; i <= opts.plotsperpage*opts.plotsperpage && it != pdfscanvaslist.end(); i++)
        {
          pagecnv->cd(i);
          (*it)->DrawClonePad();
          it++;
        }

      pgn++;
      sprintf(pgnum, "%d", pgn);
      pagecnv->Print((opts.outdir + "plots_" + pgnum + ".eps").c_str());
    }
  for (vector <TCanvas*>::iterator it = pdfscanvasratiolist.begin(); it != pdfscanvasratiolist.end();)
    {
      char numb[25];
      sprintf(numb, "ratio_%ld", it - pdfscanvasratiolist.begin());
      TCanvas * pagecnv = new TCanvas(numb, "", 0, 0, opts.resolution *opts.plotsperpage, opts.resolution * opts.plotsperpage);
      pagecnv->Divide(opts.plotsperpage, opts.plotsperpage);
      for (int i = 1; i <= opts.plotsperpage*opts.plotsperpage && it != pdfscanvasratiolist.end(); i++)
        {
          pagecnv->cd(i);
          (*it)->DrawClonePad();
          it++;
        }
      pgn++;
      snprintf(pgnum, sizeof(pgnum), "%d", pgn);
      pagecnv->Print((opts.outdir + "plots_" + pgnum + ".eps").c_str());
    }

  if (opts.twopanels || opts.threepanels)
    gStyle->SetPaperSize(opts.pagewidth, opts.pagewidth);
  else
    gStyle->SetPaperSize(opts.pagewidth*float(opts.plotsperpage)/2., opts.pagewidth*float(opts.plotsperpage)/2.);
  it = datapullscanvaslist.begin();
  for (it = datapullscanvaslist.begin(); it != datapullscanvaslist.end();)
    {
      char numb[25];
      sprintf(numb, "data_%ld", it - datapullscanvaslist.begin());
      TCanvas * pagecnv;
      if (opts.twopanels || opts.threepanels)
        {
          pagecnv = new TCanvas(numb, "", 0, 0, opts.resolution * 2, opts.resolution * 2);
          pagecnv->Divide(1, 2);
          for (int i = 1; i <= 2; i++)
            if (it != datapullscanvaslist.end())
              {
                pagecnv->cd(i);
                (*it)->DrawClonePad();
                it++;
              }
        }
      else
        {
          pagecnv = new TCanvas(numb, "", 0, 0, opts.resolution * 2, opts.resolution * 2);
          pagecnv->Divide(opts.plotsperpage, opts.plotsperpage);
          for (int i = 1; i <= opts.plotsperpage*opts.plotsperpage; i++)
            if (it != datapullscanvaslist.end())
              {
                pagecnv->cd(i);
                (*it)->DrawClonePad();
                it++;
              }
        }
      pgn++;
      snprintf(pgnum, sizeof(pgnum), "%d", pgn);
      pagecnv->Print((opts.outdir + "plots_" + pgnum + ".eps").c_str());
    }

  gStyle->SetPaperSize(opts.pagewidth, opts.pagewidth);
  it = shiftcanvaslist.begin();
  for (it = shiftcanvaslist.begin(); it != shiftcanvaslist.end(); it++)
    {
      char numb[25];
      sprintf(numb, "shift_%ld", it - shiftcanvaslist.begin());
      TCanvas * pagecnv = new TCanvas(numb, "", 0, 0, 2 * opts.resolution, (*it)->GetWindowHeight());
      (*it)->DrawClonePad();
      pgn++;
      snprintf(pgnum, sizeof(pgnum), "%d", pgn);
      pagecnv->Print((opts.outdir + "plots_" + pgnum + ".eps").c_str());
    }

  gStyle->SetPaperSize(opts.pagewidth, opts.pagewidth);
  it = chi2scancanvaslist.begin();
  for (it = chi2scancanvaslist.begin(); it != chi2scancanvaslist.end();)
    {
      char numb[25];
      sprintf(numb, "chi2scan_%ld", it - chi2scancanvaslist.begin());
      TCanvas * pagecnv;
      pagecnv = new TCanvas(numb, "", 0, 0, opts.resolution * 2, opts.resolution * 2);
      pagecnv->Divide(2, 2);
      for (int i = 1; i <= 4; i++)
        if (it != chi2scancanvaslist.end())
          {
            pagecnv->cd(i);
            (*it)->DrawClonePad();
            it++;
          }
      pgn++;
      snprintf(pgnum, sizeof(pgnum), "%d", pgn);
      pagecnv->Print((opts.outdir + "plots_" + pgnum + ".eps").c_str());
    }

  gStyle->SetPaperSize(opts.pagewidth, opts.pagewidth);
  it = chi2scangausscanvaslist.begin();
  for (it = chi2scangausscanvaslist.begin(); it != chi2scangausscanvaslist.end();)
    {
      char numb[25];
      snprintf(numb, sizeof(numb), "chi2scan_gauss_%ld", it - chi2scangausscanvaslist.begin());
      TCanvas * pagecnv;
      pagecnv = new TCanvas(numb, "", 0, 0, opts.resolution * 2, opts.resolution * 2);
      pagecnv->Divide(2, 2);
      for (int i = 1; i <= 4; i++)
	if (it != chi2scangausscanvaslist.end())
	  {
	    pagecnv->cd(i);
	    (*it)->DrawClonePad();
	    it++;
	  }
      pgn++;
      snprintf(pgnum, sizeof(pgnum), "%d", pgn);
      pagecnv->Print((opts.outdir + "plots_" + pgnum + ".eps").c_str());
    }
  
  if (opts.splitplots)
    {
      string ext = opts.ext;
      if (opts.ext == "pdf")
        ext = "eps";

      gStyle->SetPaperSize(opts.pagewidth / 2., opts.pagewidth / 2.);
      for (vector <TCanvas*>::iterator it = pdfscanvaslist.begin(); it != pdfscanvaslist.end(); it++)
        (*it)->Print((opts.outdir + (*it)->GetName() + "." + ext).c_str());
      for (vector <TCanvas*>::iterator it = pdfscanvasratiolist.begin(); it != pdfscanvasratiolist.end(); it++)
        (*it)->Print((opts.outdir + (*it)->GetName() + "." + ext).c_str());
      if (opts.twopanels || opts.threepanels)
        gStyle->SetPaperSize(opts.pagewidth, opts.pagewidth / 2.);
      for (vector <TCanvas*>::iterator it = datapullscanvaslist.begin(); it != datapullscanvaslist.end(); it++)
        (*it)->Print((opts.outdir + (*it)->GetName() + "." + ext).c_str());
      gStyle->SetPaperSize(opts.pagewidth, opts.pagewidth);
      for (vector <TCanvas*>::iterator it = shiftcanvaslist.begin(); it != shiftcanvaslist.end(); it++)
        (*it)->Print((opts.outdir + (*it)->GetName() + "." + ext).c_str());
      for (vector <TCanvas*>::iterator it = chi2scancanvaslist.begin(); it != chi2scancanvaslist.end(); it++)
        (*it)->Print((opts.outdir + (*it)->GetName() + "." + ext).c_str());

      if (opts.ext == "pdf")
        {
          for (vector <TCanvas*>::iterator it = pdfscanvaslist.begin(); it != pdfscanvaslist.end(); it++)
            system(((string)"ps2pdf -dEPSCrop " + opts.outdir + (*it)->GetName() + ".eps " + opts.outdir + (*it)->GetName() + ".pdf").c_str());
          for (vector <TCanvas*>::iterator it = pdfscanvasratiolist.begin(); it != pdfscanvasratiolist.end(); it++)
            system(((string)"ps2pdf -dEPSCrop " + opts.outdir + (*it)->GetName() + ".eps " + opts.outdir + (*it)->GetName() + ".pdf").c_str());
          for (vector <TCanvas*>::iterator it = datapullscanvaslist.begin(); it != datapullscanvaslist.end(); it++)
            system(((string)"ps2pdf -dEPSCrop " + opts.outdir + (*it)->GetName() + ".eps " + opts.outdir + (*it)->GetName() + ".pdf").c_str());
          for (vector <TCanvas*>::iterator it = shiftcanvaslist.begin(); it != shiftcanvaslist.end(); it++)
            system(((string)"ps2pdf -dEPSCrop " + opts.outdir + (*it)->GetName() + ".eps " + opts.outdir + (*it)->GetName() + ".pdf").c_str());
          for (vector <TCanvas*>::iterator it = chi2scancanvaslist.begin(); it != chi2scancanvaslist.end(); it++)
            system(((string)"ps2pdf -dEPSCrop " + opts.outdir + (*it)->GetName() + ".eps " + opts.outdir + (*it)->GetName() + ".pdf").c_str());
        }
      cout << "Multiple " << opts.ext << " plots saved in: " << (opts.outdir + "*." + opts.ext) << endl;
    }

  //Save TCanvas and TGraphs
  if (opts.root)
    {
      TFile * f = new TFile((opts.outdir + "plots.root").c_str(), "recreate");
      f->mkdir("Canvas");
      f->cd("Canvas");
      for (vector <TCanvas*>::iterator it = pdfscanvaslist.begin(); it != pdfscanvaslist.end(); it++)
        (*it)->Write();
      for (vector <TCanvas*>::iterator it = pdfscanvasratiolist.begin(); it != pdfscanvasratiolist.end(); it++)
        (*it)->Write();
      for (vector <TCanvas*>::iterator it = datapullscanvaslist.begin(); it != datapullscanvaslist.end(); it++)
        (*it)->Write();
      for (vector <TCanvas*>::iterator it = shiftcanvaslist.begin(); it != shiftcanvaslist.end(); it++)
        (*it)->Write();
      for (vector <TCanvas*>::iterator it = chi2scancanvaslist.begin(); it != chi2scancanvaslist.end(); it++)
        (*it)->Write();

      f->cd("");
      f->mkdir("Graphs");
      f->cd("Graphs");
      for (vector <TGraphAsymmErrors*>::iterator git = allgraphs.begin(); git != allgraphs.end(); git++)
        (*git)->Write();

      f->Close();
      cout << "TCanvas saved in: " << (opts.outdir + "plots.root") << endl;
    }

  //make plots file
  string format = opts.format;

  string inputfiles = "";
  for (int n = 1; n <= pgn; n++)
    {
      snprintf(pgnum, sizeof(pgnum), "%d", n);
      inputfiles = inputfiles + " " + opts.outdir + "plots_" + pgnum + ".eps";
    }

  //  if (!chi2tab)
  if (!opts.notables)
    inputfiles = inputfiles + " " + opts.outdir + "chi2.pdf";
  if (!partab)
    inputfiles = inputfiles + " " + opts.outdir + "par.pdf";
  if (!uncsummary)
    inputfiles = inputfiles + " " + opts.outdir + "unc_summary.pdf";

  //A4 is /PageSize [842 595]
  string gscommand = "gs -dBATCH -q -sDEVICE=" + format + "write -sOutputFile=" + opts.outdir + "plots." + format
    + " -dNOPAUSE -dEPSFitPage -c \"<< /PageSize [595 595] >> setpagedevice\"  -f "
    + inputfiles;

  bool makeplots = system(gscommand.c_str());
  if (!makeplots)
    cout << "Plots saved in: " << (opts.outdir + "plots." + format) << endl;
  else
    {
      cout << "ghostcript error in making: " << (opts.outdir + "plots." + format) << endl;
      cout << gscommand << endl;
    }

  //cleanup pages
  if (!opts.splitplots)
    if (pgn > 0)
      system(("rm " + inputfiles).c_str());

  return 0;
}
