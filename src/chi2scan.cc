#include "xfitter_cpp.h"

#include <string>
#include <iomanip>

//return error if LHAPDF is not enabled
#if !defined LHAPDF_ENABLED
void chi2_scan_()
{
  string msg = "S: Call to chi2_scan but LHAPDF is not enabled. Run ./configure --enable-lhapdf and link the executable";
  hf_errlog_(14060204, msg.c_str(), msg.size());
}
#elif !defined ROOT_ENABLED
void chi2_scan_()
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

void fitchi2_and_store(map <double, double> chi2, double& min, double& deltap, double& deltam, double& chi2min, string name)
{
  TGraph *chi2graph = new TGraph(chi2.size());
  int i = 0;
  for (map<double, double>::iterator it = chi2.begin(); it != chi2.end(); it++, i++)
    chi2graph->SetPoint(i, it->first, it->second);

  //2nd order fit
  TF1 *cf = new TF1("ParFit", "pol2");
  chi2graph->Fit(cf, "WQ", "", chi2graph->GetX()[0], chi2graph->GetX()[chi2graph->GetN()-1]);
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
  double min2 = parfit->GetParameter(2);
  double delta2 = parfit->GetParameter(1);
  double chi2min2 = parfit->GetParameter(0);

  //3rd order fit
  TF1 *parfit3 = new TF1("ParFit3", "pol3");
  chi2graph->Fit(parfit3, "WQ", "", chi2graph->GetX()[0], chi2graph->GetX()[chi2graph->GetN()-1]);
  double a3 = parfit3->GetParameter(3);
  double b3 = parfit3->GetParameter(2);
  double c3 = parfit3->GetParameter(1);
  double d3 = parfit3->GetParameter(0);
  double min3 = parfit3->GetMinimumX(chi2graph->GetX()[0], chi2graph->GetX()[chi2graph->GetN()-1]);
  double chi2min3 = parfit3->Eval(min3);
  TF1 *fs = new TF1("fs", "abs([3]*x**3 +[2]*x**2 + [1]*x + [0])");
  fs->SetParameter(3,a3);
  fs->SetParameter(2,b3);
  fs->SetParameter(1,c3);
  fs->SetParameter(0,d3-(chi2min3+1.));
  double deltap3 = fs->GetMinimumX(min3,chi2graph->GetX()[chi2graph->GetN()-1])-min3;
  double deltam3 = min3-fs->GetMinimumX(chi2graph->GetX()[0],min3);

  //4th order fit
  TF1 *parfit4 = new TF1("ParFit4", "pol4");
  chi2graph->Fit(parfit4, "WQ", "", chi2graph->GetX()[0], chi2graph->GetX()[chi2graph->GetN()-1]);
  double a4 = parfit4->GetParameter(4);
  double b4 = parfit4->GetParameter(3);
  double c4 = parfit4->GetParameter(2);
  double d4 = parfit4->GetParameter(1);
  double e4 = parfit4->GetParameter(0);
  double min4 = parfit4->GetMinimumX(chi2graph->GetX()[0], chi2graph->GetX()[chi2graph->GetN()-1]);
  double chi2min4 = parfit4->Eval(min4);
  TF1 *fs4 = new TF1("fs4", "abs([4]*x**4 + [3]*x**3 +[2]*x**2 + [1]*x + [0])");
  fs4->SetParameter(4,parfit4->GetParameter(4));
  fs4->SetParameter(3,parfit4->GetParameter(3));
  fs4->SetParameter(2,parfit4->GetParameter(2));
  fs4->SetParameter(1,parfit4->GetParameter(1));
  fs4->SetParameter(0,parfit4->GetParameter(0)-(chi2min4+1.));
  double deltap4 = fs4->GetMinimumX(min4,chi2graph->GetX()[chi2graph->GetN()-1])-min4;
  double deltam4 = min4-fs4->GetMinimumX(chi2graph->GetX()[0],min4);

  //5th order fit
  TF1 *parfit5 = new TF1("ParFit5", "pol5");
  chi2graph->Fit(parfit5, "WQ", "", chi2graph->GetX()[0], chi2graph->GetX()[chi2graph->GetN()-1]);
  double min5 = parfit5->GetMinimumX(chi2graph->GetX()[0], chi2graph->GetX()[chi2graph->GetN()-1]);
  double chi2min5 = parfit5->Eval(min5);
  TF1 *fs5 = new TF1("fs5", "abs([5]*pow(x,5) + [4]*x**4 + [3]*x**3 +[2]*x**2 + [1]*x + [0])");
  fs5->SetParameter(5,parfit5->GetParameter(5));
  fs5->SetParameter(4,parfit5->GetParameter(4));
  fs5->SetParameter(3,parfit5->GetParameter(3));
  fs5->SetParameter(2,parfit5->GetParameter(2));
  fs5->SetParameter(1,parfit5->GetParameter(1));
  fs5->SetParameter(0,parfit5->GetParameter(0)-(chi2min5+1.));
  double deltap5 = fs5->GetMinimumX(min5,chi2graph->GetX()[chi2graph->GetN()-1])-min5;
  double deltam5 = min5-fs5->GetMinimumX(chi2graph->GetX()[0],min5);
  
  /*
  cout << "2th order" <<endl;
  cout << "min " << min2 << "+-" << delta2 << endl;
  cout << "chi2 " << chi2min2 << endl;
  cout << "3th order" <<endl;
  cout << "min " << min3 << "+" << deltap3 << " -" << deltam3 << endl;
  cout << "chi2 " << chi2min3 << endl;
  cout << "4th order" <<endl;
  cout << "min " << min4 << "+" << deltap4 << " -" << deltam4 << endl;
  cout << "chi2 " << chi2min4 << endl;
  cout << "5th order" <<endl;
  cout << "min " << min5 << "+" << deltap5 << " -" << deltam5 << endl;
  cout << "chi2 " << chi2min5 << endl;
  */

  /*
  min = min3;
  chi2min = chi2min3;
  deltap = deltap3;
  deltam = deltam3;
  */
  min = min4;
  chi2min = chi2min4;
  deltap = deltap4;
  deltam = deltam4;
  
  //store chi2
  string outdir = string(coutdirname_.outdirname_, 128);
  outdir.erase(outdir.find_last_not_of(" ")+1, string::npos);

  string label = string(chi2scan_.label_, 64);
  if (label.find_first_of(" ") != 0)
    label.erase(label.find_last_not_of(" ")+1, string::npos);

  ofstream fchi2((outdir + "/" + name).c_str());

  fchi2 << setprecision(16);
  fchi2 << label << endl;
  fchi2 << min2 << "\t" << delta2 << "\t" << chi2min2 << endl;
  fchi2 << "ax^2+bx+c" << "\t" << a << "\t" << b << "\t" << c << endl;
  fchi2 << min3 << "\t" << deltap3 << "\t" << deltam3 << "\t" << chi2min3 << endl;
  fchi2 << "ax^3+bx^2+cx+d" << "\t" << a3 << "\t" << b3 << "\t" << c3 << "\t" << d3 << endl;
  fchi2 << min4 << "\t" << deltap4 << "\t" << deltam4 << "\t" << chi2min4 << endl;
  fchi2 << "ax^4+bx^3+cx^2+dx+e" << "\t" << a4 << "\t" << b4 << "\t" << c4 << "\t" << d4 << "\t" << e4 << endl;

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

  //check if the parameter label correspond to a parameter in xfitter
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
	      char cid[10];
	      sprintf(cid, "%d", *dit);
	      char vl[10];
	      sprintf(cid, "%f.3", *vit);
	      char tm[10];
	      sprintf(cid, "%s", *tit);
	      string msg = (string)"S: Error in chi2scan namelist: source not found for value, "  + vl + ", dataset " + cid + ", term " + tm;
	      hf_errlog_(16012001, msg.c_str(), msg.size());
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
      {
	tTEmap::iterator it = gTEmap.find(*dit);
	if (it == gTEmap.end())
	  {
	    char cid[10];
	    sprintf(cid, "%d", *dit);
	    string msg = (string)"S: Error in chi2scan namelist: dataset with ID " + cid + " has no theory expression";
	    hf_errlog_(16011801, msg.c_str(), msg.size());
	  }
	centralsources[*dit][*tit] = gTEmap[*dit]->GetTheorySource(*tit);
      }

  //Set reference PDF set if specified, and initialise theory as pseudo data if requested in MCErrors namelist
  if (lhapdfref.size() != 0)
    {
      LHAPDF::initPDFSet(lhapdfref.c_str());
      LHAPDF::initPDF(0);
      c_alphas_.alphas_ = LHAPDF::alphasPDF(boson_masses_.mz_);
    }
  //  mc_method_();
  chi2data_theory_(1);

  //variables which stores results
  double min, deltap, deltam, chi2min;
  double central, eplus, eminus;

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
	  fitchi2_and_store (chi2, min, deltap, deltam, chi2min, "chi2scan.txt");
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
	  int ModPDFErr = 0;
	  int ParPDFErr = 0;
	  getpdfunctype_heraf_(MonteCarloPDFErr, AsymHessPDFErr, SymmHessPDFErr, lhapdfset.c_str(), lhapdfset.size());
      
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

		      //In VAR PDF set determine if it is a model or parametrisation variation
		      if (pdfset == 1)
			{
			  MonteCarloPDFErr = 0;
			  AsymHessPDFErr = 0;
			  SymmHessPDFErr = 0;
			  ParPDFErr = 0;
			  ModPDFErr = 0;
			  if (iset > (nsets - chi2scan_.chi2nparvar_))
			    ParPDFErr = true;
			  else
			    ModPDFErr = true;
			}

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
	  fitchi2_and_store (chi2, min, deltap, deltam, chi2min,"chi2scan.txt");
	  xi.push_back(min);

	  //Loop on PDF sets
	  for (int iset = 1; iset <= totset; iset++)
	    {
	      double min_i, deltap_i, deltam_i, chi2min_i;
	      char chi2name[100];
	      sprintf(chi2name, "chi2scan_%d.txt", iset);
	      fitchi2_and_store (pdfchi2[iset], min_i, deltap_i, deltam_i, chi2min_i, chi2name);
	      //Need to set Scale 68
	      //min_i = min + (min_i-min)/1.64;
	      xi.push_back(min_i);
	    }

	  //Calculate error
	  central = xi[0];
	  //asymmetric hessian
	  if (AsymHessPDFErr)
	    ahessdeltaasym(xi, eplus, eminus);

	  //symmetric hessian
	  if (SymmHessPDFErr)
	    eplus = eminus = shessdelta(xi);

	  //MC errors
	  if (MonteCarloPDFErr)
	    {
	      central = median(xi);
	      deltaasym(xi, central, eplus, eminus, cl(1));
	    }
      
	  //model and parameterisation
	  if (sets == 2)
	    vardeltaasym(xi, chi2scan_.chi2nparvar_, eplus, eminus); //need to set npar

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

      cout << endl;
      cout << "Results of the chi2 scan: " << endl;
      cout << label << " = " << min << " +" << deltap << " -" << deltam << endl;
      cout << "Chi2 at minimum: " << chi2min << endl;
      if (lhapdferror)
	cout << "PDF uncertainties: " <<  central << "+" << eplus << "-" << eminus << endl;
      cout << endl;
    }

  //Fit mode
}
#endif
