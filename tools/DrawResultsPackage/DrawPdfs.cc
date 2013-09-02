#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <TError.h>
#include <TCanvas.h>
#include <TFile.h>

#include <Output.h>
#include <PdfsPainter.h>
#include <CommandParser.h>
#include <DataPainter.h>

using std::cout;
using std::cerr;
using std::endl;

using namespace std;

int main(int argc, char **argv) 
{
  //parse command line arguments
  opts = CommandParser(argc, argv);

  gErrorIgnoreLevel=1001;

  //read output directory
  vector <Output*> info_output;
  for (vector<string>::iterator it = opts.dirs.begin(); it != opts.dirs.end(); it++)
    info_output.push_back(new Output((*it).c_str()));

  //--------------------------------------------------
  //Pdf plots
  typedef map <int, vector<gstruct> > pdfmap;

  map <float, pdfmap> q2list;
  vector<string>::iterator itn = opts.labels.begin();
  for (unsigned int o = 0; o < info_output.size(); o++)
    {
      //      info_output[o]->Prepare(true);
      info_output[o]->Prepare(opts.dobands);
      info_output[o]->PreparePdf(opts.dobands);
      //loop on Q2 bins
      for (int nq2 = 0; nq2 < (info_output[o]->GetNQ2Files() / 2); nq2++) //why do I get double number of NQ2 files??!!
	{
	  pdfmap pmap;
	  if (q2list.size() > nq2)
	    pmap = q2list[info_output[o]->GetQ2Value(nq2)];

	  //loop on pdf types
	  for (unsigned int ipdf = 0; ipdf < pdflabels.size(); ipdf++)
	    {
	      //	      TGraphAsymmErrors* pdf = info_output[o]->GetPdf((Output::pdf)ipdf, nq2);
	      gstruct gs;
	      gs.graph = info_output[o]->GetPdf((Output::pdf)ipdf, nq2);
	      gs.label = (*itn);
	      pmap[ipdf].push_back(gs);
	    }
	  q2list[info_output[o]->GetQ2Value(nq2)] = pmap;
	}
      itn++;
    }

  vector <TCanvas*> pdfscanvaslist, pdfscanvasratiolist;
  for (map <float, pdfmap>::iterator qit = q2list.begin(); qit != q2list.end(); qit++)
    for (pdfmap::iterator pdfit = (*qit).second.begin(); pdfit != (*qit).second.end(); pdfit++)
      {
	pdfscanvaslist.push_back(PdfsPainter((*qit).first, (*pdfit).first, (*pdfit).second));
	pdfscanvasratiolist.push_back(PdfsRatioPainter((*qit).first, (*pdfit).first, (*pdfit).second));
      }


  //--------------------------------------------------
  //Data pulls plots
  map <int, vector<dataseth> > datamap;
  itn = opts.labels.begin(); //vector<string>::iterator 
  for (unsigned int o = 0; o < info_output.size(); o++)
    {
      for (unsigned int d = 0; d < info_output[o]->GetNsets(); d++)
	{
	  if (info_output[o]->GetSet(d)->GetNSubPlots() == 1)
	    {
	      info_output[o]->GetSet(d)->GetHistogram(0, false);
	      int id = info_output[o]->GetSet(d)->GetSetId();
	      //check bins sanity
	      vector <float> b1 = info_output[o]->GetSet(d)->getbins1(0);
	      vector <float> b2 = info_output[o]->GetSet(d)->getbins2(0);
	      vector<float>::iterator it1 = b1.begin();
	      vector<float>::iterator it2 = b2.begin();
	      bool skip = false;
	      for (; (it1+1) != b1.end(); it1++, it2++)
		if (*(it1+1) != *it2 || *it1 >= *(it1+1))
		  skip = true;
	      if (skip)
		{
		    cout << "bin inconsistency for dataset: " << info_output[o]->GetSet(d)->GetName() << endl;
		    cout << "Cannot plot data pulls, skipping" << endl;
		    continue;
		}

	      dataseth dt = dataseth(info_output[o]->GetSet(d)->GetName(),
				     info_output[o]->GetName()->Data(),
				     (*itn),
				     info_output[o]->GetSet(d)->getbins1(0),
				     info_output[o]->GetSet(d)->getbins2(0),
				     info_output[o]->GetSet(d)->getdata(0),
				     info_output[o]->GetSet(d)->getuncor(0),
				     info_output[o]->GetSet(d)->gettoterr(0),
				     info_output[o]->GetSet(d)->gettheory(0),
				     info_output[o]->GetSet(d)->gettheoryshifted(0),
				     info_output[o]->GetSet(d)->getpulls(0),
				     info_output[o]->GetSet(d)->GetXlog(0),
				     info_output[o]->GetSet(d)->GetYlog(0),
				     info_output[o]->GetSet(d)->GetXmin(0),
				     info_output[o]->GetSet(d)->GetXmax(0),
				     info_output[o]->GetSet(d)->GetXTitle(0),
				     info_output[o]->GetSet(d)->GetYTitle(0));
	      datamap[id].push_back(dt);
	    }
	}
      itn++;
    }

  vector <TCanvas*> datapullscanvaslist;
  for (map <int, vector <dataseth> >::iterator it = datamap.begin(); it != datamap.end(); it++)
    datapullscanvaslist.push_back(DataPainter((*it).first, (*it).second));



  //Save plots
  system(((string)"mkdir -p " + opts.outdir).c_str());
  vector <TCanvas*>::iterator  it = pdfscanvaslist.begin();
  (*it)->Print((opts.outdir + "Plots.eps[").c_str());
  for (vector <TCanvas*>::iterator it = pdfscanvaslist.begin(); it != pdfscanvaslist.end(); it++)
    {
      (*it)->Print((opts.outdir + "Plots.eps").c_str());
      if (opts.splitplots)
	(*it)->Print((opts.outdir + (*it)->GetName() + ".eps").c_str());
    }
  for (vector <TCanvas*>::iterator it = pdfscanvasratiolist.begin(); it != pdfscanvasratiolist.end(); it++)
    {
      (*it)->Print((opts.outdir + "Plots.eps").c_str());
      if (opts.splitplots)
	(*it)->Print((opts.outdir + (*it)->GetName() + ".eps").c_str());
    }
  for (vector <TCanvas*>::iterator it = datapullscanvaslist.begin(); it != datapullscanvaslist.end(); it++)
    {
      (*it)->Print((opts.outdir + "Plots.eps").c_str());
      if (opts.splitplots)
	(*it)->Print((opts.outdir + (*it)->GetName() + ".eps").c_str());
    }
  it = pdfscanvaslist.begin();
  (*it)->Print((opts.outdir + "Plots.eps]").c_str());

  if (opts.pdf)
    {
      cout << "Converting to pdf format..." << endl;
      system(((string)"ps2pdf " + opts.outdir + "Plots.eps " + opts.outdir + "Plots.pdf").c_str());
    }

  //Save all plots in a root file
  TFile * f = new TFile((opts.outdir + "plots.root").c_str(), "recreate");
  for (vector <TCanvas*>::iterator it = pdfscanvaslist.begin(); it != pdfscanvaslist.end(); it++)
    (*it)->Write();
  for (vector <TCanvas*>::iterator it = pdfscanvasratiolist.begin(); it != pdfscanvasratiolist.end(); it++)
    (*it)->Write();
  for (vector <TCanvas*>::iterator it = datapullscanvaslist.begin(); it != datapullscanvaslist.end(); it++)
    (*it)->Write();
  f->Close();

  return 0;
}
