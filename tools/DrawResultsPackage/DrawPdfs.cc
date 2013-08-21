#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <TError.h>
#include <TCanvas.h>

#include <Output.h>
#include <PdfsPainter.h>

using std::cout;
using std::cerr;
using std::endl;

using namespace std;

//useful trick to declare a vector "statically"
string pdflab[] = {"g", "u", "d", "uv", "dv", "#bar{U}", "#bar{D}", "Sea", "S", "C", "B"};
vector <string> pdflabels(pdflab, pdflab + sizeof(pdflab) / sizeof(string));
//be carefull! pdf names must be mapped to pdf definition in Output.h
// enum pdf{kNULL=-1, kGluon=0, kU=1, kD=2, kUv=3, kDv=4, kUb=5, kDb=6, kSea=7, kS=8, kC=9, kB=10};


int main(int argc, char **argv) 
{

  string OutputPath("output/");
  
  gErrorIgnoreLevel=1001;

  vector <string> dirs, labels;

  if (argc == 1) 
    {
      //------------------------
      cout << "program usage:\n DrawPdfs dir1[:label1] [dir2:[label2]] [...]" << endl;
      cout << " dir1 is set as reference for the ratio plots" << endl;
      //------------------------
      return -1;
    }
  
  for (int iar = 1; iar < argc; iar++)
    dirs.push_back(argv[iar]);

  if (dirs.size() > 6)
    {
      //------------------------
      cout << "program usage:\n DrawPdfs dir1[:label1] [dir2[:label2]] [...]" << endl;
      cout << " dir1 is set as reference for the ratio plots" << endl;
      cout << " maximum number of directory is 6" << endl;
      //------------------------
      return -1;
    }

  //parse dirs for lables
  for (vector<string>::iterator it = dirs.begin(); it != dirs.end(); it++)
    if ((*it).find(":") != string::npos)
      {
      labels.push_back((*it).substr((*it).find(":")+1, (*it).size() - (*it).find(":") - 1));
      (*it).erase((*it).find(":"), (*it).size());
      }
    else
      labels.push_back((*it));

  //initialize datasets
  if (dirs.size() == 1){
    OutputPath = dirs[0];
    cout << "plots are stored in: " << (dirs[0]).c_str() << endl;
  }
  else{
    OutputPath = "pdfplots/";
    cout << "plots are stored in: pdfplots/" << endl;
  }

  vector <Output*> info_output;
  for (vector<string>::iterator it = dirs.begin(); it != dirs.end(); it++)
    info_output.push_back(new Output((*it).c_str()));

  //  typedef map <string, vector<TGraphAsymmErrors*> > pdfmap;
  typedef map <string, vector<gstruct> > pdfmap;

  map <double, pdfmap> q2list;
  vector<string>::iterator itn = labels.begin();
  for (unsigned int o = 0; o < info_output.size(); o++)
    {
      info_output[o]->Prepare(true);
      info_output[o]->PreparePdf(true);
      //loop on Q2 bins
      for (int nq2 = 0; nq2 < (info_output[o]->GetNQ2Files() / 2); nq2++) //why do I get double number of NQ2 files??!!
	{
	  pdfmap pmap;
	  if (q2list.size() > nq2)
	    pmap = q2list[info_output[o]->GetQ2Value(nq2)];

	  //loop ond pdf types
	  for (unsigned int ipdf = 0; ipdf < pdflabels.size(); ipdf++)
	    {
	      //	      TGraphAsymmErrors* pdf = info_output[o]->GetPdf((Output::pdf)ipdf, nq2);
	      gstruct gs;
	      gs.graph = info_output[o]->GetPdf((Output::pdf)ipdf, nq2);
	      gs.label = (*itn);
	      pmap[pdflabels[ipdf]].push_back(gs);
	    }
	  q2list[info_output[o]->GetQ2Value(nq2)] = pmap;
	}
      itn++;
    }

  vector <TCanvas*> pdfscanvaslist, pdfscanvasratiolist;
  for (map <double, pdfmap>::iterator qit = q2list.begin(); qit != q2list.end(); qit++)
    for (pdfmap::iterator pdfit = (*qit).second.begin(); pdfit != (*qit).second.end(); pdfit++)
      {
	pdfscanvaslist.push_back(PdfsPainter((*qit).first, (*pdfit).first, (*pdfit).second));
	pdfscanvasratiolist.push_back(PdfsRatioPainter((*qit).first, (*pdfit).first, (*pdfit).second));
      }

  system(((string)"mkdir -p " + OutputPath).c_str());
  for (vector <TCanvas*>::iterator it = pdfscanvaslist.begin(); it != pdfscanvaslist.end(); it++)
    {
       string out = OutputPath + (*it)->GetName() + ".eps";
       //       string out = OutputPath + (*it)->GetName() + ".pdf";
      (*it)->Print(out.c_str());
    }

  for (vector <TCanvas*>::iterator it = pdfscanvasratiolist.begin(); it != pdfscanvasratiolist.end(); it++)
    {
       string out = OutputPath + (*it)->GetName() + ".eps";
       //       string out = OutputPath + (*it)->GetName() + ".pdf";
      (*it)->Print(out.c_str());
    }

  return 0;
}

