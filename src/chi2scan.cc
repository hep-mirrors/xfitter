#include "herafitter_cpp.h"

//return error if LHAPDF is not enabled
#if !defined LHAPDF_ENABLED
void chi2_scan_()
{
  string msg = "S: Call to chi2_scan but LHAPDF is not enabled. Run ./configure --enable-lhapdf and link the executable";
  hf_errlog_(14060204, msg.c_str(), msg.size());
}
#elif !defined ROOT_ENABLED
void chi2_scan()
{
  string msg = "S: Call to chi2_scan but ROOT library are not linked. Run ./configure with root available in your PATH";
  hf_errlog_(14062501, msg.c_str(), msg.size());
}
#elif !defined APPLGRID_ENABLED
void chi2_scan_()
{
  string msg = "S: Call to chi2_scan but ROOT library are not linked. Run ./configure with root available in your PATH";
  hf_errlog_(14062501, msg.c_str(), msg.size());
}
#else

#include <LHAPDF/LHAPDF.h>

#include "TheorEval.h"
#include "pdferrors.h"

#include <TGraph.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TF1.h>

#include <iostream>
#include <fstream>
#include <iomanip>

void fitchi2(map <double, double> chi2, double& min, double& delta, double& chi2min)
{
  TGraph *chi2graph = new TGraph(chi2.size());
  int i = 0;
  for (map<double, double>::iterator it = chi2.begin(); it != chi2.end(); it++, i++)
    chi2graph->SetPoint(i, it->first, it->second);

  TF1 *cf = new TF1("ParFit", "pol2");
  chi2graph->Fit(cf, "WQ0", "", chi2graph->GetX()[0], chi2graph->GetX()[chi2graph->GetN()-1]);
  double a = cf->GetParameter(2);
  double b = cf->GetParameter(1);
  double c = cf->GetParameter(0);
  double xc = 0;
  double sigma = 1;
  if (a > 0)
    {
      xc = -b / (2*a);
      sigma = 1. / sqrt(a);
    }
  TF1 *parfit = new TF1("ParFit", "[0]+(x-[2])**2/[1]**2");
  parfit->SetParameter(0,c);
  parfit->SetParameter(1,sigma);
  parfit->SetParameter(2,xc);
  chi2graph->Fit(parfit, "WQ", "", chi2graph->GetX()[0], chi2graph->GetX()[chi2graph->GetN()-1]);
  chi2graph->SetLineWidth(1);
  parfit->SetLineStyle(2);
  parfit->SetLineWidth(1);
  min = parfit->GetParameter(2);
  delta = parfit->GetParameter(1);
  chi2min = parfit->GetParameter(0);
}

void storechi2(map <double, double> chi2, double& min, double& delta, double& chi2min, string name)
{
  string outdir = string(coutdirname_.outdirname_, 128);
  outdir.erase(outdir.find_last_not_of(" ")+1, string::npos);

  string label = string(chi2scan_.label_, 64);
  if (label.find_first_of(" ") != 0)
    label.erase(label.find_last_not_of(" ")+1, string::npos);

  ofstream fchi2((outdir + "/" + name).c_str());

  fchi2 << label << endl;
  fchi2 << min << "\t" << delta << "\t" << chi2min << endl;

  for (map<double, double>::iterator it = chi2.begin(); it != chi2.end(); it++)
    fchi2 << it->first << "\t" << it->second << endl;

  fchi2.close();
}

void chi2_scan_()
{
  cout << endl << endl << endl;
  cout << "  -----------------------" << endl;
  cout << "  Start chi2 scan"         << endl;
  cout << "  -----------------------" << endl;
  cout << endl << endl << endl;

  //Read steering info from chi2scan fortran namelist
  string label = string(chi2scan_.label_, 64);
  if (label.find_first_of(" ") != 0)
    label.erase(label.find_last_not_of(" ")+1, string::npos);

  vector <double> values;
  for (int i = 0; i < 100; i++)
    if (chi2scan_.values_[i] == 0 && (i < 99 && chi2scan_.values_[i] == chi2scan_.values_[i+1])) //stops if two consecutive 0 are found
      break;
    else
      values.push_back(chi2scan_.values_[i]);

  //check if the parameter label correspond to a parameter in herafitter
  //else load applgrid or table scan
  vector <int> dataid;
  for (int i = 0; i < 150; i++)
    if (chi2scan_.dataid_[i] == 0)
      break;
    else
      dataid.push_back(chi2scan_.dataid_[i]);

  //Store expression vectors of terms in map: terms[dataid]
  map <int, vector <string> > terms;
  for (vector<int>::iterator dit = dataid.begin(); dit != dataid.end(); dit++)
    for (int i = 0; i < 16; i++)
      {
	string term = string(chi2scan_.term_[dit-dataid.begin()][i], 8);
	term.erase(term.find(" "), string::npos);
	if (term.size() == 0)
	  break;
	terms[*dit].push_back(term);
      }

  //Store theory sources in map: sources[values][dataid][term]
  map <double, map <int, map <string, string > > > sources;
  for (vector <double>::iterator vit = values.begin(); vit != values.end(); vit++)
    for (vector<int>::iterator dit = dataid.begin(); dit != dataid.end(); dit++)
      for (vector<string>::iterator tit = terms[*dit].begin(); tit != terms[*dit].end(); tit++)
	{
	  string source = string (chi2scan_.theorysources_[dit-dataid.begin()][tit-terms[*dit].begin()][vit-values.begin()], 1000);
	  if (source.find(" ") == string::npos)
	    {
	      cout << "Error, mismatch between number of values and number of sources in chi2scan namelist" << endl;
	      return;
	    }
	  source.erase(source.find_first_of(" "), string::npos);
	  sources[*vit][*dit][*tit] = source;
	}

  /*
  //check
  vector< vector <string> >::iterator it = sources.begin();
  for (vector <double>::iterator vit = values.begin(); vit != values.end(); vit++, it++)
    {
      vector<string>::iterator sit = (*it).begin();
      for (vector<string>::iterator tit = terms.begin(); tit != terms.end(); tit++, sit++)
	cout << *tit << " " << *sit << "  " << *vit << endl;
    }
  //endcheck
  */

  //check if lhapdf sets are specified
  string lhapdfref = string(chi2scan_.chi2lhapdfref_, 128);
  lhapdfref.erase(lhapdfref.find_first_of(" "), string::npos);
  string lhapdfset = string(chi2scan_.chi2lhapdfset_, 128);
  lhapdfset.erase(lhapdfset.find_first_of(" "), string::npos);
  string lhapdfvarset = string(chi2scan_.chi2lhapdfvarset_, 128);
  lhapdfvarset.erase(lhapdfvarset.find_first_of(" "), string::npos);
  bool lhapdferror = chi2scan_.pdferrors_;

  string outdir = string(coutdirname_.outdirname_, 128);
  outdir.erase(outdir.find_last_not_of(" ")+1, string::npos);

  map <int, map <string, string> > centralsources;
  for (vector<int>::iterator dit = dataid.begin(); dit != dataid.end(); dit++)
    for (vector<string>::iterator tit = terms[*dit].begin(); tit != terms[*dit].end(); tit++)
      centralsources[*dit][*tit] = gTEmap[*dit]->GetTheorySource(*tit);

  //Set reference PDF set if specified, and initialise theory as pseudo data if requested in MCErrors namelist
  if (lhapdfref.size() != 0)
    {
      LHAPDF::initPDFSet(lhapdfref.c_str());
      LHAPDF::initPDF(0);
      c_alphas_.alphas_ = LHAPDF::alphasPDF(boson_masses_.mz_);
    }
  //  mc_method_();
  chi2data_theory_(1);

  if (lhapdfset.size() != 0)
    {
      //simple scan mode, to estimate experimental uncertainty
      if (!lhapdferror)
	{
	  LHAPDF::initPDFSet(lhapdfset.c_str());
	  LHAPDF::initPDF(0);
	  map <double, double> chi2; //map of parameters value and chi2 values

	  //Loop on parameters points
	  for (vector <double>::iterator vit = values.begin(); vit != values.end(); vit++)
	    {
	      for (vector<int>::iterator dit = dataid.begin(); dit != dataid.end(); dit++)
		for (vector<string>::iterator tit = terms[*dit].begin(); tit != terms[*dit].end(); tit++)
		  gTEmap[*dit]->ChangeTheorySource(*tit, sources[*vit][*dit][*tit]);
	      double chi2tot = chi2data_theory_(2);
	      char chi2c[20];
	      sprintf(chi2c, "%.2f", chi2tot);
	      chi2[*vit] = chi2tot;
	    }
	  double min, delta, chi2min;
	  fitchi2 (chi2, min, delta, chi2min);
	  cout << min << "  " << delta << endl;
	  storechi2 (chi2, min, delta, chi2min, "chi2scan.txt");
	}
      else       //lhapdferror mode
	{
	  int totset = 0; //total number of PDF members
	  int sets = 1;
	  if (lhapdfvarset != "")
	    sets = 2;
	  int MonteCarloPDFErr = 0;
	  int AsymHessPDFErr = 0;
	  int SymmHessPDFErr = 0;
	  getpdfunctype_heraf_(lhapdfset.c_str(), MonteCarloPDFErr, AsymHessPDFErr, SymmHessPDFErr);
      
	  //Build the chi2 map
	  map <double, double> chi2; //central PDF map of parameters value and chi2 values
	  map <int, map <double, double> > pdfchi2; //map of PDF members and map of parameters value and chi2
	  
	  //Loop on parameters points
	  for (vector <double>::iterator vit = values.begin(); vit != values.end(); vit++)
	    {
	      for (vector<int>::iterator dit = dataid.begin(); dit != dataid.end(); dit++)
		for (vector<string>::iterator tit = terms[*dit].begin(); tit != terms[*dit].end(); tit++)
		  gTEmap[*dit]->ChangeTheorySource(*tit, sources[*vit][*dit][*tit]);
	      
	      int cset = 1; //counter on pdf variations;
	      for (int pdfset = 0; pdfset < sets; pdfset++)
		{
		  if (pdfset ==  0 && sets > 1)
		    LHAPDF::initPDFSet(lhapdfset.c_str());
		  else if (pdfset ==  1)
		    {
		      if (lhapdfvarset == "")
			continue;
		      LHAPDF::initPDFSet(lhapdfvarset.c_str());
		    }
		  //Number of PDF members
		  int nsets = LHAPDF::numberPDF();
		  cout << "Number of PDF members for this set: " << nsets << endl;
		  if (vit == values.begin())
		    totset += nsets;

		  //Loop on PDF sets
		  for (int iset = 0; iset <= nsets; iset++)
		    {
		      //Init PDF member and alphas(MZ)
		      LHAPDF::initPDF(iset);
		      c_alphas_.alphas_ = LHAPDF::alphasPDF(boson_masses_.mz_);
		      double chi2tot = chi2data_theory_(2);
		      char chi2c[20];
		      sprintf(chi2c, "%.2f", chi2tot);
		      //		      cout << setw(15) << "chi2=" << chi2c 
		      //			   << setw(15) << "ndf=" << cfcn_.ndfmini_ 
		      //			   << endl;
		      if (iset == 0)
			chi2[*vit] = chi2tot;
		      else
			{
			  pdfchi2[cset][*vit] = chi2tot;
			  cset++;
			}
		    }
		}
	    }
      
	  vector <double> xi;
	  //Central PDF
	  double min0, delta0, chi2min0;
	  fitchi2 (chi2, min0, delta0, chi2min0);
	  cout << min0 << "  " << delta0 << endl;
	  storechi2 (chi2, min0, delta0, chi2min0, "chi2scan.txt");
	  xi.push_back(min0);

	  //Loop on PDF sets
	  for (int iset = 1; iset <= totset; iset++)
	    {
	      double min, delta, chi2min;
	      fitchi2 (pdfchi2[iset], min, delta, chi2min);
	      //Need to set Scale 68
	      //min = min0 + (min-min0)/1.64;
	      char chi2name[100];
	      sprintf(chi2name, "chi2scan_%d.txt", iset);
	      storechi2 (pdfchi2[iset], min, delta, chi2min, chi2name);
	      xi.push_back(min);
	    }

	  //Calculate error
	  double central = xi[0];
	  double eplus, eminus;
	  //asymmetric hessian
	  if (AsymHessPDFErr)
	    ahessdeltaasym(xi, eplus, eminus);

	  //symmetric hessian
	  if (SymmHessPDFErr)
	    eplus = eminus = shessdelta(xi);

	  //MC errors
	  if (MonteCarloPDFErr)
	    {
	      central = mean(xi);
	      deltaasym(xi, central, eplus, eminus, cl(1));
	    }
      
	  //model and parameterrisation
	  if (sets == 2)
	    vardeltaasym(xi, 0, eplus, eminus); //need to set npar

	  cout << "Result: " <<  central << "+" << eplus << "-" << eminus << endl;
	}
      //print fittedresults.txt with nominal grid and test PDF
      for (vector<int>::iterator dit = dataid.begin(); dit != dataid.end(); dit++)
	for (vector<string>::iterator tit = terms[*dit].begin(); tit != terms[*dit].end(); tit++)
	  gTEmap[*dit]->ChangeTheorySource(*tit, centralsources[*dit][*tit]);
      LHAPDF::initPDFSet(lhapdfset.c_str());
      LHAPDF::initPDF(0);
      string fname = outdir + "/Results.txt";
      fopen_(85, fname.c_str(), fname.size());
      double chi2tot = chi2data_theory_(3);
      fclose_(85);
    }

  //Fit mode
}
#endif
