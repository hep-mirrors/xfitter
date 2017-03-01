#include "xfitter_cpp.h"
#include "dimensions.h"

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

struct point
{
  double thc;                //central theory
  vector <double> th_asym_p; //asymmetric plus variations
  vector <double> th_asym_m; //asymmetric minus variations
  vector <double> th_hess_s; //hessian symmetric variations
  vector <double> th_par;    //parametrisation variations
  double th_env_p;           //up parametrisation envelope
  double th_env_m;           //down parametrisation envelope
  vector <double> th_mc;     //monte carlo variations
  double th_mc_mean;         //monte carlo mean
  double th_mc_var;          //monte carlo var
  vector <double> th_scale_p; //scale plus variations
  vector <double> th_scale_m; //scale minus variations
  double th_err_up;          //theory error up
  double th_err_dn;          //theory error down
};


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
  string outdir = string(coutdirname_.outdirname_, 256);
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
  //read the label assigned to the scan parameter
  string label = string(chi2scan_.label_, 64);
  if (label.find_first_of(" ") != 0)
    label.erase(label.find_last_not_of(" ")+1, string::npos);

  //read the parameters values
  vector <double> values;
  for (int i = 0; i < NCHI2POINTS_C; i++)
    if (chi2scan_.values_[i] == 0 && ((i+1) < NCHI2POINTS_C && chi2scan_.values_[i] == chi2scan_.values_[i+1])) //stops if two consecutive 0 are found
      break;
    else
      values.push_back(chi2scan_.values_[i]);

  //To implement: check if the parameter label correspond to a parameter in xfitter, else load applgrid or table scan

  //read the data ids
  vector <int> dataid;
  for (int i = 0; i < NSET_C; i++)
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
	      sprintf(vl, "%.3f", *vit);
	      string msg = (string)"S: Error in chi2scan namelist: source not found for value, "  + vl + ", dataset " + cid + ", term " + (*tit);
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
  bool lhapdfprofile = chi2scan_.pdfprofile_;
  bool scaleprofile = chi2scan_.scaleprofile_;
  //  bool scaleprofile = true;

  string outdir = string(coutdirname_.outdirname_, 256);
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

  if (lhapdfset.size() == 0)
    {
      //Fit mode: perform a PDF fit at each value of the parameter (not implemented)
      string msg = "S: Call to chi2_scan but lhapdfset is not specified. Fit mode is not yet implemented";
      hf_errlog_(16022601, msg.c_str(), msg.size());
    }

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

	  //double chi2tot = chi2data_theory_(2);
	  char vl[10];
	  sprintf(vl, "%.3f", *vit);
	  string fname = outdir + "/Results_" + vl + ".txt";
	  fopen_(85, fname.c_str(), fname.size());
	  double chi2tot = chi2data_theory_(3);
	  fclose_(85);
	  writefittedpoints_(); //write out fittedresults.txt file
	  bool mv = system(((string)"mv " + outdir + "/fittedresults.txt " 
			    + outdir + "/fittedresults_" + vl + ".txt").c_str());
	  
	  char chi2c[20];
	  sprintf(chi2c, "%.1f", chi2tot);
	  cout << setw(15) << (label + "=") << *vit
	       << setw(15) << "chi2=" << chi2c 
	       << setw(15) << "ndf=" << cfcn_.ndfmini_ 
	       << endl;
	  chi2[*vit] = chi2tot;
	}
      fitchi2_and_store (chi2, min, deltap, deltam, chi2min, "chi2scan.txt");
    }
  else if(lhapdferror && !lhapdfprofile) //lhapdferror mode
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
	  
      //Loop on parameters values
      for (vector <double>::iterator vit = values.begin(); vit != values.end(); vit++)
	{
	  //change Theory source to the current parameter value
	  for (vector<int>::iterator dit = dataid.begin(); dit != dataid.end(); dit++)
	    for (vector<string>::iterator tit = terms[*dit].begin(); tit != terms[*dit].end(); tit++)
	      gTEmap[*dit]->ChangeTheorySource(*tit, sources[*vit][*dit][*tit]);
	      
	  //Start loop on PDF sets and members
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

	      //Loop on PDF members
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

		  double chi2tot;
		  if (pdfset == 0 && iset == 0)
		    {
		      //char chi2c[20];
		      //sprintf(chi2c, "%.1f", chi2tot);
		      //cout << setw(15) << (label + "=") << *vit
		      //<< setw(15) << "chi2=" << chi2c 
		      //<< setw(15) << "ndf=" << cfcn_.ndfmini_ 
		      //<< endl;
		      char vl[10];
		      sprintf(vl, "%.3f", *vit);
		      string fname = outdir + "/Results_" + vl + ".txt";
		      fopen_(85, fname.c_str(), fname.size());
		      chi2tot = chi2data_theory_(3);
		      fclose_(85);
		      writefittedpoints_(); //write out fittedresults.txt file
		      bool mv = system(((string)"mv " + outdir + "/fittedresults.txt " 
					+ outdir + "/fittedresults_" + vl + ".txt").c_str());
		    }
		  else
		    chi2tot = chi2data_theory_(2);
		  
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
  else if(lhapdfprofile) //lhapdferror profile mode
    {
      int npoints = cndatapoints_.npoints_;

      int nsysexp = systema_.nsys_; //initial number of systematic uncertainties, before addition of PDF uncertainties
      
      map <int, point> pointsmap;
      map <double, double> chi2; //map of parameters value and chi2 values
      
      int totset = 0; //total number of PDF members
      int sets = 1;
      if (lhapdfvarset != "")
	sets = 2;
      
      //Loop on parameters values
      for (vector <double>::iterator vit = values.begin(); vit != values.end(); vit++)
	{
	  int nsysloc = nsysexp; //initial number of systematic uncertainties, before addition of PDF uncertainties

	  //Reset number of systematic uncertainties
	  systema_.nsys_ = nsysexp;
	  systematicsflags_.resetcommonsyst_ = true;
	  
	  int MonteCarloPDFErr = 0;
	  int AsymHessPDFErr = 0;
	  int SymmHessPDFErr = 0;
	  int ModPDFErr = 0;
	  int ParPDFErr = 0;
	  getpdfunctype_heraf_(MonteCarloPDFErr, AsymHessPDFErr, SymmHessPDFErr, lhapdfset.c_str(), lhapdfset.size());
	  
	  //change Theory source to the current parameter value
	  for (vector<int>::iterator dit = dataid.begin(); dit != dataid.end(); dit++)
	    for (vector<string>::iterator tit = terms[*dit].begin(); tit != terms[*dit].end(); tit++)
	      gTEmap[*dit]->ChangeTheorySource(*tit, sources[*vit][*dit][*tit]);
	  
	  //clean up points
	  for (map <int, point>::iterator  pit = pointsmap.begin(); pit != pointsmap.end(); pit++)
	    {
	      (*pit).second.th_asym_p.clear();
	      (*pit).second.th_asym_m.clear();
	      (*pit).second.th_hess_s.clear();
	      (*pit).second.th_par.clear();   
	      (*pit).second.th_mc.clear();    
	      (*pit).second.th_scale_p.clear();
	      (*pit).second.th_scale_m.clear();
	    }

	  //Start loop on PDF sets and members
	  int cset = 0; //counter on pdf variations;
	  int isys = 1; //counter on PDF error index;  (not used)
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

	      //Loop on PDF members
	      for (int iset = 0; iset <= nsets; iset++)
		{
		  //skip central member of the VAR set
		  if (pdfset == 1 && iset == 0)
		    continue;
		  
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
		  
		  //Store all PDF variations for each data point
		  for (int i = 0; i < npoints; i++)
		    if (cset == 0) //central PDF
		      pointsmap[i].thc = c_theo_.theo_[i];
		    else if (iset != 0) // PDF variations
		      {
			if (MonteCarloPDFErr)
			  pointsmap[i].th_mc.push_back(c_theo_.theo_[i]);
			else if (SymmHessPDFErr)
			  pointsmap[i].th_hess_s.push_back(c_theo_.theo_[i]);
			else if (AsymHessPDFErr)
			  if ((cset%2) == 1)
			    pointsmap[i].th_asym_p.push_back(c_theo_.theo_[i]);
			  else
			    pointsmap[i].th_asym_m.push_back(c_theo_.theo_[i]);
			else if (ModPDFErr)
			  if ((cset%2) == 1)
			    pointsmap[i].th_asym_p.push_back(c_theo_.theo_[i]);
			  else
			    pointsmap[i].th_asym_m.push_back(c_theo_.theo_[i]);
			else if (ParPDFErr)
			  pointsmap[i].th_par.push_back(c_theo_.theo_[i]);
		      }
	  
		  if (iset != 0)
		    if (MonteCarloPDFErr || SymmHessPDFErr || ParPDFErr)
		      isys++;
		    else if ((iset%2) == 0) //set the same index for Up and Down variation of asymmetric PDF errors
		      isys++;
		  cset++;

		} //end loop on PDF members
	    } //end loop on PDF sets

	  //In QCD scales profiling mode, add two nuisance parameters for the renormalisation and factorisation scales
	  if (scaleprofile)
	    {
	      LHAPDF::initPDFSet(lhapdfset.c_str());
	      LHAPDF::initPDF(0);
	      c_alphas_.alphas_ = LHAPDF::alphasPDF(boson_masses_.mz_);
	      map <int, int> iordmap;
	      map <int, double> murmap;
	      map <int, double> mufmap;
	      for (vector<int>::iterator dit = dataid.begin(); dit != dataid.end(); dit++)
		{
		  int iord;
		  double mur0;
		  double muf0;
		  gTEmap[*dit]->GetOrdScales(iord, mur0, muf0);
		  iordmap[*dit] = iord;
		  murmap[*dit] = mur0;
		  mufmap[*dit] = muf0;
		}

	      double factor = 2.;
	      double chi2tot;

	      //mur*2
	      for (vector<int>::iterator dit = dataid.begin(); dit != dataid.end(); dit++)
		gTEmap[*dit]->SetOrdScales(iordmap[*dit], murmap[*dit]*factor, mufmap[*dit]);
	      chi2tot = chi2data_theory_(2);
	      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
		pointsmap[i].th_scale_p.push_back(c_theo_.theo_[i]);
	      //mur*0.5
	      for (vector<int>::iterator dit = dataid.begin(); dit != dataid.end(); dit++)
		gTEmap[*dit]->SetOrdScales(iordmap[*dit], murmap[*dit]/factor, mufmap[*dit]);
	      chi2tot = chi2data_theory_(2);
	      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
		pointsmap[i].th_scale_m.push_back(c_theo_.theo_[i]);

	      //muf*2
	      for (vector<int>::iterator dit = dataid.begin(); dit != dataid.end(); dit++)
		gTEmap[*dit]->SetOrdScales(iordmap[*dit], murmap[*dit], mufmap[*dit]*factor);
	      chi2tot = chi2data_theory_(2);
	      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
		pointsmap[i].th_scale_p.push_back(c_theo_.theo_[i]);
	      //muf*0.5
	      for (vector<int>::iterator dit = dataid.begin(); dit != dataid.end(); dit++)
		gTEmap[*dit]->SetOrdScales(iordmap[*dit], murmap[*dit], mufmap[*dit]/factor);
	      chi2tot = chi2data_theory_(2);
	      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
		pointsmap[i].th_scale_m.push_back(c_theo_.theo_[i]);

	      //restore nominal scale
	      for (vector<int>::iterator dit = dataid.begin(); dit != dataid.end(); dit++)
		gTEmap[*dit]->SetOrdScales(iordmap[*dit], murmap[*dit], mufmap[*dit]);
	    }
	  
	  //Add PDF uncertainties as nuisance parameters to the chi2
	  
	  //Monte Carlo replica PDF uncertainties
	  int totmc = pointsmap.begin()->second.th_mc.size();;
	  if (totmc > 0)
	    {

	      //MC replica mean and rms
	      for (map <int, point>::iterator  pit = pointsmap.begin(); pit != pointsmap.end(); pit++)
		{
		  vector <double> xi;
		  for (vector <double>::iterator mcit = pit->second.th_mc.begin(); mcit != pit->second.th_mc.end(); mcit++)
		    xi.push_back(*mcit);
		  pit->second.th_mc_mean = mean(xi);
		  pit->second.th_mc_var = rms(xi);
		}
	      char num[10];
	      sprintf (num, "%d", totmc);
	      string msg = (string) "I: Found " + num + " Monte Carlo PDF uncertainties variations";
	      hf_errlog_(25051401, msg.c_str(), msg.size());

	      //Evaluate MC covariance matrix
	      unsigned int dim = npoints;
	      double *covmx = new  double [dim*dim];
	      for (map <int, point>::iterator  pit1 = pointsmap.begin(); pit1 != pointsmap.end(); pit1++)
		for (map <int, point>::iterator  pit2 = pointsmap.begin(); pit2 != pointsmap.end(); pit2++)
		  {
		    int i = pit1->first;
		    int j = pit2->first;
		    if (pit1->second.th_mc.size() != pit2->second.th_mc.size())
		      {
			string msg = (string)"S: Error: inconsistent number of MC replica per point";
			hf_errlog_(25051402, msg.c_str(), msg.size());
		      }

		    covmx[j*dim+i] = 0;


		    for (int mc = 0; mc < totmc; mc++)
		      covmx[j*dim+i] += (pit1->second.th_mc[mc] - pit1->second.th_mc_mean) * (pit2->second.th_mc[mc] - pit2->second.th_mc_mean) / (double)totmc;
		  }

	      //Evaluate covariance matrix to nuisance parameters conversion      
	      double *beta_from_covmx = new double [dim*npoints];
	      double alpha_from_covmx[dim];

	      int ncorr = 0;	
	      getnuisancefromcovar_(dim,npoints,npoints,
				    covmx,beta_from_covmx,0,
				    ncorr,alpha_from_covmx,0);


	      if ( ncorr + nsysloc > NSYSMAX_C ) {
		char csys[6];
		sprintf( csys,"%i", ncorr + nsysloc) ;
		string msg = (string) "S: Too many systematic sources, increase NSYSMAX_C to " + csys + " at least in  include/dimensions.h and recompile";
		hf_errlog_(15111901,msg.c_str(), msg.size());

	      }

	      for (int j = 0; j < ncorr; j++)
		{
		  for (int i = 0; i < npoints; i++)
		    {
		      systema_.beta_[i][nsysloc] = beta_from_covmx[i*dim+j]/pointsmap[i].th_mc_mean;
		      systasym_.omega_[i][nsysloc] = 0;
		      if (clhapdf_.scale68_)
			{
			  systema_.beta_[i][nsysloc] /= 1.64;
			  systasym_.omega_[i][nsysloc] /= 1.64;
			}
		    }
		  nsysloc += 1;
		}
	      delete covmx;
	      delete beta_from_covmx;
	    }
  
	  //Asymmetric PDF uncertainties, including hessian and model
	  int totasym = 0;
	  for (map <int, point>::iterator  pit = pointsmap.begin(); pit != pointsmap.end(); pit++)
	    if (pit->second.th_asym_p.size() != pit->second.th_asym_m.size())
	      {
		string msg = (string)"S: Error: inconsistent number of positive and negative asymmetric PDF variations, check your NPARVAR setting";
		hf_errlog_(25051403, msg.c_str(), msg.size());
	      }
	    else
	      totasym = pit->second.th_asym_p.size(); //include both model and parametrisation uncertainties
	  if (totasym > 0)
	    {
	      char num[10];
	      sprintf (num, "%d", totasym);
	      string msg = (string) "I: Found " + num + " asymmetric PDF uncertainties variations";
	      hf_errlog_(25051404, msg.c_str(), msg.size());
	    }
	  //Add PDF uncertainties
	  for (int j = 0; j < totasym; j++)
	    {
	      for (int i = 0; i < npoints; i++)
		{
		  //account for the sign flip due to applying a theory variation as a shift to the data
		  systasym_.betaasym_[i][0][nsysloc] = -(pointsmap[i].th_asym_p[j] - pointsmap[i].thc) / pointsmap[i].thc;
		  systasym_.betaasym_[i][1][nsysloc] = -(pointsmap[i].th_asym_m[j] - pointsmap[i].thc) / pointsmap[i].thc;
		  systema_.beta_[i][nsysloc] = 0.5*(systasym_.betaasym_[i][0][nsysloc] - systasym_.betaasym_[i][1][nsysloc]);
		  systasym_.omega_[i][nsysloc] = 0.5*(systasym_.betaasym_[i][0][nsysloc] + systasym_.betaasym_[i][1][nsysloc]);
		  
		  if (clhapdf_.scale68_)
		    {
		      systema_.beta_[i][nsysloc]  /= 1.64;
		      systasym_.betaasym_[i][1][nsysloc] /= 1.64;
		      systasym_.betaasym_[i][0][nsysloc] /= 1.64;
		      systasym_.omega_[i][nsysloc] /= 1.64;
		    }
		}
	      nsysloc += 1;
	    }

	  //Symmetric hessian uncertainties
	  int tothess_s = 0;
	  for (map <int, point>::iterator  pit = pointsmap.begin(); pit != pointsmap.end(); pit++)
	    tothess_s = pit->second.th_hess_s.size();
	  if (tothess_s > 0)
	    {
	      char num[10];
	      sprintf (num, "%d", tothess_s);
	      string msg = (string) "I: Found " + num + " Symmetric hessian PDF uncertainties variations";
	      hf_errlog_(25051405, msg.c_str(), msg.size());
	    }
	  for (int j = 0; j < tothess_s; j++)
	    {
	      for (int i = 0; i < npoints; i++)
		{
		  systema_.beta_[i][nsysloc] = (pointsmap[i].th_hess_s[j] - pointsmap[i].thc) / pointsmap[i].thc;
		  systasym_.omega_[i][nsysloc] = 0;
		  if (clhapdf_.scale68_)
		    {
		      systema_.beta_[i][nsysloc]  /= 1.64;
		      systasym_.omega_[i][nsysloc] /= 1.64;
		    }
		}
	      nsysloc += 1;
	    }

	  //Parametrisation uncertainties
	  int totpar = 0;
	  for (map <int, point>::iterator  pit = pointsmap.begin(); pit != pointsmap.end(); pit++)
	    totpar = pit->second.th_par.size();
	  if (totpar > 0)
	    {
	      char num[10];
	      sprintf (num, "%d", totpar);
	      string msg = (string) "I: Found " + num + " Parametrisation uncertainties";
	      hf_errlog_(25051406, msg.c_str(), msg.size());
	    }
	  //make envelope;
	  for (map <int, point>::iterator  pit = pointsmap.begin(); pit != pointsmap.end(); pit++)
	    {
	      vector <double> xi;
	      xi.push_back(pit->second.thc);
	      for (vector <double>::iterator it = pit->second.th_par.begin(); it != pit->second.th_par.end(); it++)
		xi.push_back(*it);

	      double pareplus, pareminus;
	      deltaenvelope(xi, pareplus, pareminus);

	      pit->second.th_env_p = pit->second.thc + pareplus;
	      pit->second.th_env_m = pit->second.thc - pareminus;
	    }
	  if (totpar > 0)
	    {
	      for (int i = 0; i < npoints; i++)
		{
		  //account for the sign flip due to applying a theory variation as a shift to the data
		  systasym_.betaasym_[i][0][nsysloc] = -(pointsmap[i].th_env_p - pointsmap[i].thc) / pointsmap[i].thc;
		  systasym_.betaasym_[i][1][nsysloc] = -(pointsmap[i].th_env_m - pointsmap[i].thc) / pointsmap[i].thc;
		  systema_.beta_[i][nsysloc] = 0.5*(systasym_.betaasym_[i][0][nsysloc] - systasym_.betaasym_[i][1][nsysloc]);
		  systasym_.omega_[i][nsysloc] = 0.5*(systasym_.betaasym_[i][0][nsysloc] + systasym_.betaasym_[i][1][nsysloc]);

		  //treat as symmetric
		  //systema_.beta_[i][nsysloc] = -(pointsmap[i].th_env_p - pointsmap[i].th_env_m) / 2./ pointsmap[i].thc;
		  //systasym_.omega_[i][nsysloc] = 0;

		  if (clhapdf_.scale68_)
		    {
		      systema_.beta_[i][nsysloc]  /= 1.64;
		      systasym_.betaasym_[i][1][nsysloc] /= 1.64;
		      systasym_.betaasym_[i][0][nsysloc] /= 1.64;
		      systasym_.omega_[i][nsysloc] /= 1.64;
		    }
		}
	      nsysloc += 1;
	    }

	  // Set PDF nuisance parameter name
	  for (int j = systema_.nsys_; j < nsysloc; j++)
	    {
	      char nuispar[64];
	      sprintf (nuispar, "PDF_nuisance_param_%02d", j+1 - systema_.nsys_);
	      int len = strlen(nuispar);
	      memset (nuispar+len,' ',64-len);
	      strcpy(systema_.system_[j],nuispar);
	    }


	  int nsyspdf = nsysloc;
	  //Asymmetric scale variations
	  int totscale = 0;
	  for (map <int, point>::iterator  pit = pointsmap.begin(); pit != pointsmap.end(); pit++)
	    if (pit->second.th_scale_p.size() != pit->second.th_scale_m.size())
	      {
		string msg = (string)"S: Error: inconsistent number of positive and negative asymmetric scale variations";
		hf_errlog_(2016030201, msg.c_str(), msg.size());
	      }
	    else
	      totscale = pit->second.th_scale_p.size();
	  if (totscale > 0)
	    {
	      char num[10];
	      sprintf (num, "%d", totscale);
	      string msg = (string) "I: Found " + num + " asymmetric scale variations";
	      hf_errlog_(2016030202, msg.c_str(), msg.size());
	    }
	  //Add scale variations
	  for (int j = 0; j < totscale; j++)
	    {
	      for (int i = 0; i < npoints; i++)
		{
		  //account for the sign flip due to applying a theory variation as a shift to the data
		  systasym_.betaasym_[i][0][nsysloc] = -(pointsmap[i].th_scale_p[j] - pointsmap[i].thc) / pointsmap[i].thc;
		  systasym_.betaasym_[i][1][nsysloc] = -(pointsmap[i].th_scale_m[j] - pointsmap[i].thc) / pointsmap[i].thc;
		  systema_.beta_[i][nsysloc] = 0.5*(systasym_.betaasym_[i][0][nsysloc] - systasym_.betaasym_[i][1][nsysloc]);
		  systasym_.omega_[i][nsysloc] = 0.5*(systasym_.betaasym_[i][0][nsysloc] + systasym_.betaasym_[i][1][nsysloc]);
		  if (clhapdf_.scale68_)
		    {
		      systema_.beta_[i][nsysloc]  /= 1.64;
		      systasym_.betaasym_[i][1][nsysloc] /= 1.64;
		      systasym_.betaasym_[i][0][nsysloc] /= 1.64;
		      systasym_.omega_[i][nsysloc] /= 1.64;
		    }
		}
	      nsysloc += 1;
	    }

	  // Set scale nuisance parameter name
	  for (int j = nsyspdf; j < nsysloc; j++)
	    {
	      char nuispar[64];
	      sprintf (nuispar, "scale_nuisance_param_%02d", j+1 - nsyspdf);
	      int len = strlen(nuispar);
	      memset (nuispar+len,' ',64-len);
	      strcpy(systema_.system_[j],nuispar);
	    }
	  
	  //Set central theory value in fortran common block
	  double theo_cent[NTOT_C];
	  for (map <int, point>::iterator  pit = pointsmap.begin(); pit != pointsmap.end(); pit++)
	    theo_cent[pit->first] = pit->second.thc;
	  bool symm = false;
	  //	  if (totmc || tothess_s)
	  //	    symm = true;
	  writetheoryfiles_(nsysloc-systema_.nsys_, theo_cent, symm);


	  LHAPDF::initPDF(0);
	  c_alphas_.alphas_ = LHAPDF::alphasPDF(boson_masses_.mz_);

	  cout << "-------------------------------------------" << endl;
	  cout << "Chi2 test on central prediction with PDF uncertainties:" << endl;

	  for (int i = systema_.nsys_; i < nsysloc; i++)
	    {
	      sysmeas_.n_syst_meas_[i] = npoints; //PDF systematic uncertainties apply to all points
	      for (int j = 0; j < npoints; j++)
		sysmeas_.syst_meas_idx_[i][j] = j + 1;
	      systscal_.sysscalingtype_[i] = 1;  //Apply linear scaling to PDF uncertainties
	      //      systscal_.sysscalingtype_[i] = 0;  //No scaling for PDF uncertainties
	      csysttype_.isysttype_[i] = 2; // THEORY
	      //systscal_.sysform_[i] = 4; //External (Minuit Fit)
	      //systscal_.sysform_[i] = 3; //Offset (Minuit Fit)
	    }
	  
	  //Add the PDF uncertainties to the total number of systematic uncertainties
	  systema_.nsys_ = nsysloc;
	  systematicsflags_.resetcommonsyst_ = true;

	  //double chi2tot = chi2data_theory_(2);
	  char vl[10];
	  sprintf(vl, "%.3f", *vit);
	  string fname = outdir + "/Results_" + vl + ".txt";
	  fopen_(85, fname.c_str(), fname.size());
	  double chi2tot = chi2data_theory_(3);
	  fclose_(85);
	  writefittedpoints_(); //write out fittedresults.txt file
	  bool mv = system(((string)"mv " + outdir + "/fittedresults.txt " 
			    + outdir + "/fittedresults_" + vl + ".txt").c_str());
	  
	  char chi2c[20];
	  sprintf(chi2c, "%.1f", chi2tot);
	  cout << setw(15) << (label + "=") << *vit
	       << setw(15) << "chi2=" << chi2c 
	       << setw(15) << "ndf=" << cfcn_.ndfmini_ 
	       << endl;
	  
	  chi2[*vit] = chi2tot;
	} //end loop on parameter values

      //Central PDF
      fitchi2_and_store (chi2, min, deltap, deltam, chi2min,"chi2scan.txt");
    }
  
  //Pick up the value closest to the minimum
  double closestval = *(values.begin());
  for (vector <double>::iterator vit = values.begin(); vit != values.end(); vit++)
    if (abs(*vit - min) < abs(closestval - min))
      closestval = *vit;

  char vl[10];
  sprintf(vl, "%.3f", closestval);
  bool cp = system(((string)"cp " + outdir + "/fittedresults_" + vl + ".txt "
		    + outdir + "/fittedresults.txt").c_str());
  cp = system(((string)"cp " + outdir + "/Results_" + vl + ".txt "
	       + outdir + "/Results.txt").c_str());

  /*
  //load theory sources corresponding to the closest value
  for (vector<int>::iterator dit = dataid.begin(); dit != dataid.end(); dit++)
    for (vector<string>::iterator tit = terms[*dit].begin(); tit != terms[*dit].end(); tit++)
      gTEmap[*dit]->ChangeTheorySource(*tit, sources[closestval][*dit][*tit]);
  */

  /*
  //load theory sources corresponding to the initial value
  for (vector<int>::iterator dit = dataid.begin(); dit != dataid.end(); dit++)
    for (vector<string>::iterator tit = terms[*dit].begin(); tit != terms[*dit].end(); tit++)
      gTEmap[*dit]->ChangeTheorySource(*tit, centralsources[*dit][*tit]);
  */

  //print fittedresults.txt and Results.txt with nominal grid and test PDF
  /*
  LHAPDF::initPDFSet(lhapdfset.c_str());
  LHAPDF::initPDF(0);
  c_alphas_.alphas_ = LHAPDF::alphasPDF(boson_masses_.mz_);
  string fname = outdir + "/Results.txt";
  fopen_(85, fname.c_str(), fname.size());
  double chi2tot = chi2data_theory_(3);
  fclose_(85);
  */


  //******************************************//
  //Store PDF members for plots
  if (lhapdferror)
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
      
      //Start loop on PDF sets and members
      int cset = 0; //counter on pdf variations;
      int isys = 1; //counter on PDF error index;
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
	  totset += nsets;
	  
	  //Loop on PDF members
	  for (int iset = 0; iset <= nsets; iset++)
	    {
	      //skip central member of the VAR set
	      if (pdfset == 1 && iset == 0)
		continue;
	      
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

	      char tag[10];
	      sprintf (tag, "_%04d", cset);
	      string filename = outdir + "/pdfs_q2val_";
	      if (iset != 0)
		{
		  if (MonteCarloPDFErr)
		    sprintf (tag, "mc%03d", isys); // Monte Carlo prefix "mc"
		  else if (SymmHessPDFErr)
		    sprintf (tag, "s%02d", isys); // Hessian prefix "s"
		  else if (AsymHessPDFErr)
		    sprintf (tag, "s%02d", isys); // Hessian prefix "s"
		  else if (ModPDFErr)
		    sprintf (tag, "m%02d", isys); // Model prefix "m"
		  else if (ParPDFErr)
		    sprintf (tag, "p%02d", isys); // Parameter prefix "p"
		  filename += tag;
		  
		  if (MonteCarloPDFErr || SymmHessPDFErr || ParPDFErr)
		    filename += "s_";   //Symmetric suffix
		  else
		    if ((cset%2) == 1)
		      filename += "p_";  //Asymmetric minus suffix
		    else
		      filename += "m_";  //Asymmetric plus suffix
		}
	      filename += " ";
	      
	      store_pdfs_(filename.c_str(), filename.size());
	      if (cset == 0)
		{
		  fill_c_common_();
		  print_lhapdf6_();
		}
	      else
		save_data_lhapdf6_(iset);

	      if (iset != 0)
		if (MonteCarloPDFErr || SymmHessPDFErr || ParPDFErr)
		  isys++;
		else if ((iset%2) == 0) //set the same index for Up and Down variation of asymmetric PDF errors
		  isys++;
	      cset++;
	    }
	}
    }
  //******************************************//

  cout << endl;
  cout << "Results of the chi2 scan: " << endl;
  cout << label << " = " << min << " +" << deltap << " -" << deltam << endl;
  cout << "Chi2 at minimum: " << chi2min << "  "  << "ndf=" << (cfcn_.ndfmini_-1) << endl;
  if (lhapdferror && ! lhapdfprofile)
    cout << "PDF uncertainties: " <<  central << "+" << eplus << "-" << eminus << endl;
  cout << endl;
}
#endif
