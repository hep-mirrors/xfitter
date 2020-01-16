#include "xfitter_cpp.h"
#include "dimensions.h"
#include "Profiler.h"
#include "xfitter_steer.h"
#include "BaseEvolution.h"

#include <string>
#include <iomanip>

#include <TError.h>

////return error if LHAPDF is not enabled
//#if !defined LHAPDF_ENABLED
//void alphas_scan_()
//{
//  string msg = "S: Call to chi2_scan but LHAPDF is not enabled. Run ./configure --enable-lhapdf and link the executable";
//  hf_errlog_(14060204, msg.c_str(), msg.size());
//}
//#elif !defined ROOT_ENABLED
//void alphas_scan_()
//{
//  string msg = "S: Call to chi2_scan but ROOT library are not linked. Run ./configure with root available in your PATH";
//  hf_errlog_(14062501, msg.c_str(), msg.size());
//}
//#elif !defined APPLGRID_ENABLED
//void alphas_scan_()
//{
//  string msg = "S: Call to chi2_scan but ROOT library are not linked. Run ./configure with root available in your PATH";
//  hf_errlog_(14062501, msg.c_str(), msg.size());
//}
//#else

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

void getpdfunctype_as(int& MonteCarloPDFErr, int& AsymHessPDFErr, int& SymmHessPDFErr, xfitter::BaseEvolution* evol)
{
  string errortype = evol->getPropertyS("ErrorType");
  if (errortype == "hessian" )
    {
      AsymHessPDFErr = true;
      SymmHessPDFErr = false;
      MonteCarloPDFErr = false;
    }
  if (errortype == "symmhessian" )
    {
      AsymHessPDFErr = false;
      SymmHessPDFErr = true;
      MonteCarloPDFErr = false;
    }
  if (errortype == "replicas" )
    {
      AsymHessPDFErr = false;
      SymmHessPDFErr = false;
      MonteCarloPDFErr = true;
    }
}

namespace asscan
{
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
    gErrorIgnoreLevel=1001;
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
    parfit3->SetParameter(3,0);
    parfit3->SetParameter(2,a);
    parfit3->SetParameter(1,b);
    parfit3->SetParameter(0,c);
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
    parfit4->SetParameter(3,a3);
    parfit4->SetParameter(2,b3);
    parfit4->SetParameter(1,c3);
    parfit4->SetParameter(0,d3);
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
    string outdir = string(coutdirname_.outdirname, 128);
    outdir.erase(outdir.find_last_not_of(" ")+1, string::npos);

    //label assigned to the scan parameter
    string label = "#alpha_{s}";

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

    //clean up memory
    delete chi2graph;
    delete cf;
    delete parfit;
    delete parfit3;
    delete parfit4;
    delete parfit5;
    delete fs;
    delete fs4;
    delete fs5;
  }


  void decompose(map <int, map <int, map <double, double> > > &systchi2, double value);
  void decompose_fits(map <int, map <int, map <double, double> > > systchi2, double min, vector <double> &deltapi, vector <double> &deltami);
}

using namespace asscan;

void alphas_scan_()
{
  cout << endl << endl << endl;
  cout << "  -----------------------" << endl;
  cout << "  Start alphas scan"         << endl;
  cout << "  -----------------------" << endl;
  cout << endl << endl << endl;

  //label assigned to the scan parameter
  string label = "#alpha_{s}";

  //Read steering info from alphasscan fortran namelist

  //check if lhapdf sets are specified
  string alphaslhapdf = string(alphasscan_.alphaslhapdf_, 128);       //PDF set for the alphas scan
  string lhapdfref = string(alphasscan_.aslhapdfref_, 128);           //PDF set for the pseudodata
  string lhapdfset = string(alphasscan_.aslhapdfset_, 128);           //PDF set for the PDF uncertainties
  string lhapdfvarset = string(alphasscan_.aslhapdfvarset_, 128);     //PDF set for the VAR PDF uncertainties

  //Remove trailing spaces from the fortran strings
  alphaslhapdf.erase(alphaslhapdf.find_first_of(" "), string::npos);
  lhapdfref.erase(lhapdfref.find_first_of(" "), string::npos);
  lhapdfset.erase(lhapdfset.find_first_of(" "), string::npos);
  lhapdfvarset.erase(lhapdfvarset.find_first_of(" "), string::npos);

  bool lhapdfprofile = alphasscan_.aspdfprofile_;
  bool scaleprofile = alphasscan_.asscaleprofile_;
  bool decomposition = true;

  //Reduce LHAPDF verbosity
  LHAPDF::Info& cfg = LHAPDF::getConfig();
  cfg.set_entry("Verbosity", 0);

  // get evolution
  auto evol=xfitter::get_evolution("proton-LHAPDF");
  YAML::Node gNode=XFITTER_PARS::getEvolutionNode("proton-LHAPDF");
  
  //Read the alphas values for the scan from the alphaslhapdf PDF set
  cout << "-------------------------------------------" << endl;
  cout << "Read the alphas values from " << alphaslhapdf << endl;
  vector <double> values;
  //LHAPDF::initPDFSet(alphaslhapdf.c_str());
  gNode["set"] = alphaslhapdf;
  gNode["member"] = 0;
  evol->atConfigurationChange();
  //int nsets = LHAPDF::numberPDF();
  int nsets = evol->getPropertyI("NumMembers")-1;
  for (int iset = 0; iset <= nsets; iset++)
    {
      //Init PDF member
      //LHAPDF::initPDF(iset);
      gNode["member"] = iset;
      evol->atConfigurationChange();
      //double as = LHAPDF::alphasPDF(boson_masses_.Mz);
      double as = evol->getAlphaS(boson_masses_.Mz);
      values.push_back(as);
    }
  cout << "Found " << values.size() << " alphas variations" << alphaslhapdf << endl;
  
  string outdir = string(coutdirname_.outdirname, 128);
  outdir.erase(outdir.find_last_not_of(" ")+1, string::npos);

  //Set reference PDF set if specified, and initialise theory as pseudo data if requested in MCErrors namelist
  if (lhapdfref.size() != 0)
    {
      //LHAPDF::initPDFSet(lhapdfref.c_str());
      //LHAPDF::initPDF(0);
      //c_alphas_.alphas = LHAPDF::alphasPDF(boson_masses_.Mz);
      //clhapdf_.ilhapdfset = 0;
      //Dyturbo::pdfname = lhapdfref;
      //Dyturbo::pdfmember = 0;
      gNode["set"] = lhapdfref;
      gNode["member"] = 0;
      evol->atConfigurationChange();
      c_alphas_.alphas = evol->getAlphaS(boson_masses_.Mz);
    }
  //  mc_method_();
  chi2data_theory_(1);

  map <int, point> pointsmap;               //This map store PDF variations for all data points
  int npoints = cndatapoints_.npoints;     //number of data points
  int nsysexp = systema_.nsys;             //initial number of systematic uncertainties, before addition of PDF and scale uncertainties
  int nsysloc = nsysexp;                              //current number of systematic uncertainties

  //Save the central theory
  if (lhapdfprofile || scaleprofile)
    {
      //Set up the central PDF set for PDF and/or scale variations
      //LHAPDF::initPDFSet(lhapdfset.c_str());
      //LHAPDF::initPDF(0);
      //c_alphas_.alphas = LHAPDF::alphasPDF(boson_masses_.Mz);
      //clhapdf_.ilhapdfset = 0;
      //Dyturbo::pdfname = lhapdfset;
      //Dyturbo::pdfmember = 0;
      gNode["set"] = lhapdfset;
      gNode["member"] = 0;
      evol->atConfigurationChange();
      c_alphas_.alphas = evol->getAlphaS(boson_masses_.Mz);

      double chi2tot = chi2data_theory_(2);

      char chi2c[500];
      sprintf(chi2c, "%.2f", chi2tot);
      cout << setw(20) << "Central PDF set: "
	   << setw(15) << "chi2/ndf = " << chi2c << "/" << cfcn_.ndfmini
	   << endl;
	      
      //Store the central theory for each data point
      for (int i = 0; i < npoints; i++)
	pointsmap[i].thc = c_theo_.theo[i];
    }
  
  //If requested, add PDF uncertainties from the uncertainty-PDF-set 
  if (lhapdfprofile)
    {
      cout << "-------------------------------------------" << endl;
      cout << "Evaluate PDF uncertainties for PDF profiling" << endl;
      
      int totset = 0; //total number of PDF members
      int sets = 1;
      if (lhapdfvarset != "")
	sets = 2;

      int MonteCarloPDFErr = 0;
      int AsymHessPDFErr = 0;
      int SymmHessPDFErr = 0;
      int ModPDFErr = 0;
      int ParPDFErr = 0;
      getpdfunctype_as(MonteCarloPDFErr, AsymHessPDFErr, SymmHessPDFErr, evol);

      //Start loop on PDF sets and members
      int cset = 1; //counter on pdf variations;
      int isys = 1; //counter on PDF error index;  (not used)
      for (int pdfset = 0; pdfset < sets; pdfset++)
	{
	  if (pdfset ==  0 && sets > 1)
	    //LHAPDF::initPDFSet(lhapdfset.c_str());
	    gNode["set"] = lhapdfset;
	  else if (pdfset ==  1)
	    {
	      if (lhapdfvarset == "")
		continue;
	      //LHAPDF::initPDFSet(lhapdfvarset.c_str());
	      //strcpy(clhapdf_.lhapdfset,lhapdfvarset.c_str());
	      //Dyturbo::pdfname = lhapdfvarset;
	      gNode["set"] = lhapdfvarset;
	    }
	  gNode["member"] = 0;
	  evol->atConfigurationChange();

	  //Number of PDF members
	  //int nsets = LHAPDF::numberPDF();
	  int nsets = evol->getPropertyI("NumMembers")-1;
	  cout << "Number of PDF members for this set: " << nsets << endl;
	  totset += nsets;

	  //Loop on PDF members
	  for (int iset = 1; iset <= nsets; iset++) //skip central member
	    {
	      //Init PDF member and alphas(MZ)
	      //LHAPDF::initPDF(iset);
	      //c_alphas_.alphas = LHAPDF::alphasPDF(boson_masses_.Mz);
	      //clhapdf_.ilhapdfset = iset;
	      //Dyturbo::pdfmember = iset;
	      gNode["member"] = iset;
	      evol->atConfigurationChange();
	      c_alphas_.alphas = evol->getAlphaS(boson_masses_.Mz);
		  
	      //In VAR PDF set determine if it is a model or parametrisation variation
	      if (pdfset == 1)
		{
		  MonteCarloPDFErr = 0;
		  AsymHessPDFErr = 0;
		  SymmHessPDFErr = 0;
		  ParPDFErr = 0;
		  ModPDFErr = 0;
		  if (iset > (nsets - chi2scan_.chi2nparvar))
		    ParPDFErr = true;
		  else
		    ModPDFErr = true;
		}

	      double chi2tot = chi2data_theory_(2);

	      char chi2c[500];
	      sprintf(chi2c, "%.2f", chi2tot);
	      cout << setw(20) << "PDF set number: " << setw(5) << iset 
		   << setw(15) << "chi2/ndf = " << chi2c << "/" << cfcn_.ndfmini
		   << endl;
	      
	      //Store all PDF variations for each data point
	      for (int i = 0; i < npoints; i++)
		{
		  if (MonteCarloPDFErr)
		    pointsmap[i].th_mc.push_back(c_theo_.theo[i]);
		  else if (SymmHessPDFErr)
		    pointsmap[i].th_hess_s.push_back(c_theo_.theo[i]);
		  else if (AsymHessPDFErr)
		    if ((cset%2) == 1)
		      pointsmap[i].th_asym_p.push_back(c_theo_.theo[i]);
		    else
		      pointsmap[i].th_asym_m.push_back(c_theo_.theo[i]);
		  else if (ModPDFErr)
		    if ((cset%2) == 1)
		      pointsmap[i].th_asym_p.push_back(c_theo_.theo[i]);
		    else
		      pointsmap[i].th_asym_m.push_back(c_theo_.theo[i]);
		  else if (ParPDFErr)
		    pointsmap[i].th_par.push_back(c_theo_.theo[i]);
		}
	  
	      if (MonteCarloPDFErr || SymmHessPDFErr || ParPDFErr)
		isys++;
	      else if ((iset%2) == 0) //set the same index for Up and Down variation of asymmetric PDF errors
		isys++;
	      
	      cset++;
	    } //end loop on PDF members
	} //end loop on PDF sets
    }
  
  //In QCD scales profiling mode, add nuisance parameters for the renormalisation, factorisation, and resummation scales
  if (scaleprofile)
    {
      //Set up the central PDF set for scale variations (same set as for PDF uncertainties)
      //LHAPDF::initPDFSet(lhapdfset.c_str());
      //LHAPDF::initPDF(0);
      //c_alphas_.alphas = LHAPDF::alphasPDF(boson_masses_.Mz);
      //clhapdf_.ilhapdfset = 0;
      //Dyturbo::pdfname = lhapdfset;
      //Dyturbo::pdfmember = 0;
      gNode["set"] = lhapdfset;
      gNode["member"] = 0;
      evol->atConfigurationChange();
      c_alphas_.alphas = evol->getAlphaS(boson_masses_.Mz);

      double factor = 2.;
      double chi2tot;

      /*
      bool mv;
      int nvar = 0;
      //mur*2
       for (map <int, TheorEval* >::iterator tit = gTEmap.begin(); tit != gTEmap.end(); tit++)
	tit->second->SetOrdScales(iordmap[tit->first], murmap[tit->first]*factor, mufmap[tit->first], muresmap[tit->first]);
      chi2tot = chi2data_theory_(2);
      char chi2c[500];
      sprintf(chi2c, "%.2f", chi2tot);
      cout << setw(20) << "mur = " << factor << ": "
	   << setw(15) << "chi2/ndf = " << chi2c << "/" << cfcn_.ndfmini
	   << endl;
      char tag[10]; sprintf (tag, "_%04d", nvar);
      writefittedpoints_(); //write out fittedresults.txt file
      mv = system(((string)"mv " + outdir + "/fittedresults.txt " + outdir + "/fittedresults.txt_scale" + tag).c_str());
      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
	pointsmap[i].th_scale_p.push_back(c_theo_.theo[i]);
      nvar++;
      //mur*0.5
      for (map <int, TheorEval* >::iterator tit = gTEmap.begin(); tit != gTEmap.end(); tit++)
	tit->second->SetOrdScales(iordmap[tit->first], murmap[tit->first]/factor, mufmap[tit->first], muresmap[tit->first]);
      chi2tot = chi2data_theory_(2);
      chi2c[500];
      sprintf(chi2c, "%.2f", chi2tot);
      cout << setw(20) << "mur = " << 1./factor << ": "
	   << setw(15) << "chi2/ndf = " << chi2c << "/" << cfcn_.ndfmini
	   << endl;
      tag[10]; sprintf (tag, "_%04d", nvar);
      writefittedpoints_(); //write out fittedresults.txt file
      mv = system(((string)"mv " + outdir + "/fittedresults.txt " + outdir + "/fittedresults.txt_scale" + tag).c_str());
      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
	pointsmap[i].th_scale_m.push_back(c_theo_.theo[i]);
      nvar++;
      
      //muf*2
      for (map <int, TheorEval* >::iterator tit = gTEmap.begin(); tit != gTEmap.end(); tit++)
	tit->second->SetOrdScales(iordmap[tit->first], murmap[tit->first], mufmap[tit->first]*factor, muresmap[tit->first]);
      chi2tot = chi2data_theory_(2);
      chi2c[500];
      sprintf(chi2c, "%.2f", chi2tot);
      cout << setw(20) << "muf = " << factor << ": "
	   << setw(15) << "chi2/ndf = " << chi2c << "/" << cfcn_.ndfmini
	   << endl;
      tag[10]; sprintf (tag, "_%04d", nvar);
      writefittedpoints_(); //write out fittedresults.txt file
      mv = system(((string)"mv " + outdir + "/fittedresults.txt " + outdir + "/fittedresults.txt_scale" + tag).c_str());
      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
	pointsmap[i].th_scale_p.push_back(c_theo_.theo[i]);
      nvar++;
      //muf*0.5
      for (map <int, TheorEval* >::iterator tit = gTEmap.begin(); tit != gTEmap.end(); tit++)
	tit->second->SetOrdScales(iordmap[tit->first], murmap[tit->first], mufmap[tit->first]/factor, muresmap[tit->first]);
      chi2tot = chi2data_theory_(2);
      chi2c[500];
      sprintf(chi2c, "%.2f", chi2tot);
      cout << setw(20) << "muf = " << 1./factor << ": "
	   << setw(15) << "chi2/ndf = " << chi2c << "/" << cfcn_.ndfmini
	   << endl;
      tag[10]; sprintf (tag, "_%04d", nvar);
      writefittedpoints_(); //write out fittedresults.txt file
      mv = system(((string)"mv " + outdir + "/fittedresults.txt " + outdir + "/fittedresults.txt_scale" + tag).c_str());
      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
	pointsmap[i].th_scale_m.push_back(c_theo_.theo[i]);
      nvar++;
      
      //mures*2
      for (map <int, TheorEval* >::iterator tit = gTEmap.begin(); tit != gTEmap.end(); tit++)
	tit->second->SetOrdScales(iordmap[tit->first], murmap[tit->first], mufmap[tit->first], muresmap[tit->first]*factor);
      chi2tot = chi2data_theory_(2);
      chi2c[500];
      sprintf(chi2c, "%.2f", chi2tot);
      cout << setw(20) << "mures = " << factor << ": "
	   << setw(15) << "chi2/ndf = " << chi2c << "/" << cfcn_.ndfmini
	   << endl;
      tag[10]; sprintf (tag, "_%04d", nvar);
      writefittedpoints_(); //write out fittedresults.txt file
      mv = system(((string)"mv " + outdir + "/fittedresults.txt " + outdir + "/fittedresults.txt_scale" + tag).c_str());
      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
	pointsmap[i].th_scale_p.push_back(c_theo_.theo[i]);
      nvar++;
      //mures*0.5
      for (map <int, TheorEval* >::iterator tit = gTEmap.begin(); tit != gTEmap.end(); tit++)
	tit->second->SetOrdScales(iordmap[tit->first], murmap[tit->first], mufmap[tit->first], muresmap[tit->first]/factor);
      chi2tot = chi2data_theory_(2);
      chi2c[500];
      sprintf(chi2c, "%.2f", chi2tot);
      cout << setw(20) << "mures = " << 1./factor << ": "
	   << setw(15) << "chi2/ndf = " << chi2c << "/" << cfcn_.ndfmini
	   << endl;
      tag[10]; sprintf (tag, "_%04d", nvar);
      writefittedpoints_(); //write out fittedresults.txt file
      mv = system(((string)"mv " + outdir + "/fittedresults.txt " + outdir + "/fittedresults.txt_scale" + tag).c_str());
      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
	pointsmap[i].th_scale_m.push_back(c_theo_.theo[i]);
      nvar++;
      
      //restore nominal scale
      for (map <int, TheorEval* >::iterator tit = gTEmap.begin(); tit != gTEmap.end(); tit++)
	tit->second->SetOrdScales(iordmap[tit->first], murmap[tit->first], mufmap[tit->first], muresmap[tit->first]);
      */

      //Consider global reaction-specific scale parameters
      double *murappl = 0;
      double *mufappl = 0;
      if (XFITTER_PARS::gParameters.find("APPLgrid/muR") != XFITTER_PARS::gParameters.end())
	murappl = XFITTER_PARS::getParamD("APPLgrid/muR");
      if (XFITTER_PARS::gParameters.find("APPLgrid/muF") != XFITTER_PARS::gParameters.end())
	mufappl = XFITTER_PARS::getParamD("APPLgrid/muF");
      double mur0appl, muf0appl;
      if (murappl)
	mur0appl = *murappl;
      if (mufappl)
	muf0appl = *mufappl;

      double *murhathor = 0;
      double *mufhathor = 0;
      if (XFITTER_PARS::gParameters.find("Hathor/muR") != XFITTER_PARS::gParameters.end())
	murhathor = XFITTER_PARS::getParamD("Hathor/muR");
      if (XFITTER_PARS::gParameters.find("Hathor/muF") != XFITTER_PARS::gParameters.end())
	mufhathor = XFITTER_PARS::getParamD("Hathor/muF");
      double mur0hathor, muf0hathor;
      if (murhathor)
	mur0hathor = *murhathor;
      if (mufhathor)
	muf0hathor = *mufhathor;

      double *murdyt = 0;
      double *mufdyt = 0;
      double *muresdyt = 0;
      if (XFITTER_PARS::gParameters.find("DYTurbo/muR") != XFITTER_PARS::gParameters.end())
	murdyt = XFITTER_PARS::getParamD("DYTurbo/muR");
      if (XFITTER_PARS::gParameters.find("DYTurbo/muF") != XFITTER_PARS::gParameters.end())
	mufdyt = XFITTER_PARS::getParamD("DYTurbo/muF");
      if (XFITTER_PARS::gParameters.find("DYTurbo/muRes") != XFITTER_PARS::gParameters.end())
	muresdyt = XFITTER_PARS::getParamD("DYTurbo/muRes");
      double mur0dyt, muf0dyt, mures0dyt;
      if (murdyt)
	mur0dyt = *murdyt;
      if (mufdyt)
	muf0dyt = *mufdyt;
      if (muresdyt)
	mures0dyt = *muresdyt;

      //mur*2
      if (murappl) *murappl = mur0appl*factor;
      if (murhathor) *murhathor = mur0hathor*factor;
      if (murdyt) *murdyt = mur0dyt*factor;
      xfitter::updateAtConfigurationChange();
      //update_theory_iteration_();
      chi2data_theory_(2);
      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
	pointsmap[i].th_scale_p.push_back(c_theo_.theo[i]);

      //mur*0.5
      if (murappl) *murappl = mur0appl/factor;
      if (murhathor) *murhathor = mur0hathor/factor;
      if (murdyt) *murdyt = mur0dyt/factor;
      xfitter::updateAtConfigurationChange();
      //update_theory_iteration_();
      chi2data_theory_(2);
      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
	pointsmap[i].th_scale_m.push_back(c_theo_.theo[i]);

      //restore nominal scale
      if (murappl) *murappl = mur0appl;
      if (murhathor) *murhathor = mur0hathor;
      if (murdyt) *murdyt = mur0dyt;
      
      //muf*2
      if (mufappl) *mufappl = muf0appl*factor;
      if (mufhathor) *mufhathor = muf0hathor*factor;
      if (mufdyt) *mufdyt = muf0dyt*factor;
      xfitter::updateAtConfigurationChange();
      //update_theory_iteration_();
      chi2data_theory_(2);
      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
	pointsmap[i].th_scale_p.push_back(c_theo_.theo[i]);
      
      //muf*0.5
      if (mufappl) *mufappl = muf0appl/factor;
      if (mufhathor) *mufhathor = muf0hathor/factor;
      if (mufdyt) *mufdyt = muf0dyt/factor;
      xfitter::updateAtConfigurationChange();
      //update_theory_iteration_();
      chi2data_theory_(2);
      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
	pointsmap[i].th_scale_m.push_back(c_theo_.theo[i]);

      //restore nominal scale
      if (mufappl) *mufappl = muf0appl;
      if (mufhathor) *mufhathor = muf0hathor;
      if (mufdyt) *mufdyt = muf0dyt;

      //mures*2
      if (muresdyt) *muresdyt = mures0dyt*factor;
      xfitter::updateAtConfigurationChange();
      //update_theory_iteration_();
      chi2data_theory_(2);
      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
	pointsmap[i].th_scale_p.push_back(c_theo_.theo[i]);
      
      //mures*0.5
      if (muresdyt) *muresdyt = mures0dyt/factor;
      xfitter::updateAtConfigurationChange();
      //update_theory_iteration_();
      chi2data_theory_(2);
      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
	pointsmap[i].th_scale_m.push_back(c_theo_.theo[i]);

      //restore nominal scale
      if (muresdyt) *muresdyt = mures0dyt;

      xfitter::updateAtConfigurationChange();
      
    }
  

  //Add PDF uncertainties as nuisance parameters to the chi2
  if (lhapdfprofile)
    {
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
		  systema_.beta[i][nsysloc] = beta_from_covmx[i*dim+j]/pointsmap[i].th_mc_mean;
		  systasym_.omega[i][nsysloc] = 0;
		  if (clhapdf_.scale68)
		    {
		      systema_.beta[i][nsysloc] /= 1.64;
		      systasym_.omega[i][nsysloc] /= 1.64;
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
	  sysmeas_.n_syst_meas[nsysloc] = npoints;
	  for (int i = 0; i < npoints; i++)
	    {
	      sysmeas_.syst_meas_idx[nsysloc][i] = i + 1;

	      //account for the sign flip due to applying a theory variation as a shift to the data
	      systasym_.betaasym[i][0][nsysloc] = -(pointsmap[i].th_asym_p[j] - pointsmap[i].thc) / pointsmap[i].thc;
	      systasym_.betaasym[i][1][nsysloc] = -(pointsmap[i].th_asym_m[j] - pointsmap[i].thc) / pointsmap[i].thc;
	      systema_.beta[i][nsysloc] = 0.5*(systasym_.betaasym[i][0][nsysloc] - systasym_.betaasym[i][1][nsysloc]);
	      systasym_.omega[i][nsysloc] = 0.5*(systasym_.betaasym[i][0][nsysloc] + systasym_.betaasym[i][1][nsysloc]);
		  
	      systasym_.lasymsyst[nsysloc] = true;
	      
	      if (clhapdf_.scale68)
		{
		  systema_.beta[i][nsysloc]  /= 1.64;
		  systasym_.betaasym[i][1][nsysloc] /= 1.64;
		  systasym_.betaasym[i][0][nsysloc] /= 1.64;
		  systasym_.omega[i][nsysloc] /= 1.64;
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
	  sysmeas_.n_syst_meas[nsysloc] = npoints;
	  for (int i = 0; i < npoints; i++)
	    {
	      sysmeas_.syst_meas_idx[nsysloc][i] = i + 1;

	      systema_.beta[i][nsysloc] = (pointsmap[i].th_hess_s[j] - pointsmap[i].thc) / pointsmap[i].thc;
	      systasym_.omega[i][nsysloc] = 0;

	      systasym_.lasymsyst[nsysloc] = false;

	      if (clhapdf_.scale68)
		{
		  systema_.beta[i][nsysloc]  /= 1.64;
		  systasym_.omega[i][nsysloc] /= 1.64;
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
	  sysmeas_.n_syst_meas[nsysloc] = npoints;
	  for (int i = 0; i < npoints; i++)
	    {
	      sysmeas_.syst_meas_idx[nsysloc][i] = i + 1;

	      //account for the sign flip due to applying a theory variation as a shift to the data
	      systasym_.betaasym[i][0][nsysloc] = -(pointsmap[i].th_env_p - pointsmap[i].thc) / pointsmap[i].thc;
	      systasym_.betaasym[i][1][nsysloc] = -(pointsmap[i].th_env_m - pointsmap[i].thc) / pointsmap[i].thc;
	      systema_.beta[i][nsysloc] = 0.5*(systasym_.betaasym[i][0][nsysloc] - systasym_.betaasym[i][1][nsysloc]);
	      systasym_.omega[i][nsysloc] = 0.5*(systasym_.betaasym[i][0][nsysloc] + systasym_.betaasym[i][1][nsysloc]);

	      systasym_.lasymsyst[nsysloc] = true;

	      //treat as symmetric
	      //systema_.beta[i][nsysloc] = -(pointsmap[i].th_env_p - pointsmap[i].th_env_m) / 2./ pointsmap[i].thc;
	      //systasym_.omega[i][nsysloc] = 0;

	      if (clhapdf_.scale68)
		{
		  systema_.beta[i][nsysloc]  /= 1.64;
		  systasym_.betaasym[i][1][nsysloc] /= 1.64;
		  systasym_.betaasym[i][0][nsysloc] /= 1.64;
		  systasym_.omega[i][nsysloc] /= 1.64;
		}
	    }
	  nsysloc += 1;
	}

      // Set PDF nuisance parameter name
      for (int j = systema_.nsys; j < nsysloc; j++)
	{
	  char nuispar[64];
	  sprintf (nuispar, "PDF_nuisance_param_%02d", j+1 - systema_.nsys);
	  int len = strlen(nuispar);
	  memset (nuispar+len,' ',64-len);
	  strcpy(systema_.system[j],nuispar);
	}
    }

  int nsyspdf = nsysloc;
  if (scaleprofile)
    {
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
	  sysmeas_.n_syst_meas[nsysloc] = npoints;
	  for (int i = 0; i < npoints; i++)
	    {
	      sysmeas_.syst_meas_idx[nsysloc][i] = i + 1;

	      //account for the sign flip due to applying a theory variation as a shift to the data
	      systasym_.betaasym[i][0][nsysloc] = -(pointsmap[i].th_scale_p[j] - pointsmap[i].thc) / pointsmap[i].thc;
	      systasym_.betaasym[i][1][nsysloc] = -(pointsmap[i].th_scale_m[j] - pointsmap[i].thc) / pointsmap[i].thc;
	      systema_.beta[i][nsysloc] = 0.5*(systasym_.betaasym[i][0][nsysloc] - systasym_.betaasym[i][1][nsysloc]);
	      systasym_.omega[i][nsysloc] = 0.5*(systasym_.betaasym[i][0][nsysloc] + systasym_.betaasym[i][1][nsysloc]);

	      systasym_.lasymsyst[nsysloc] = true;

	      if (clhapdf_.scale68)
		{
		  systema_.beta[i][nsysloc]  /= 1.64;
		  systasym_.betaasym[i][1][nsysloc] /= 1.64;
		  systasym_.betaasym[i][0][nsysloc] /= 1.64;
		  systasym_.omega[i][nsysloc] /= 1.64;
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
	  strcpy(systema_.system[j],nuispar);
	}
    }

  if (lhapdfprofile || scaleprofile)
    {
      //Set central theory value in fortran common block
      double theo_cent[NTOT_C];
      for (map <int, point>::iterator  pit = pointsmap.begin(); pit != pointsmap.end(); pit++)
	theo_cent[pit->first] = pit->second.thc;
      bool symm = false;
      //	  if (totmc || tothess_s)
      //	    symm = true;
      writetheoryfiles_(nsysloc-systema_.nsys, theo_cent, symm);

      cout << "-------------------------------------------" << endl;
      cout << "Add PDF and scale uncertainties and compute chi-square test on the central prediction:" << endl;

      //LHAPDF::initPDF(0);
      //c_alphas_.alphas = LHAPDF::alphasPDF(boson_masses_.Mz);
      //clhapdf_.ilhapdfset = 0;
      //Dyturbo::pdfmember = 0;
      gNode["member"] = 0;
      evol->atConfigurationChange();
      c_alphas_.alphas = evol->getAlphaS(boson_masses_.Mz);

      for (int i = systema_.nsys; i < nsysloc; i++)
	{
	  sysmeas_.n_syst_meas[i] = npoints; //PDF systematic uncertainties apply to all points
	  for (int j = 0; j < npoints; j++)
	    sysmeas_.syst_meas_idx[i][j] = j + 1;
	  systscal_.sysscalingtype[i] = 1;  //Apply linear scaling to PDF uncertainties
	  //systscal_.sysscalingtype[i] = 0;  //No scaling for PDF uncertainties
	  csysttype_.isysttype[i] = 2; // THEORY
	  //systscal_.sysform[i] = 4; //External (Minuit Fit)
	  //systscal_.sysform[i] = 3; //Offset (Minuit Fit)
	  if (((nsyspdf-i) <= clhapdf_.nremovepriors) && ((nsyspdf-i) > 0))
	    {
	      char nuispar[64];
	      sprintf (nuispar, "PDF_nuisance_param_%02d", i+1 - systema_.nsys);
	      cout << "Remove prior for syst " << nuispar << endl;
	      string msg = (string) "I: Remove prior for systematic " + nuispar;
	      hf_errlog_(15082401, msg.c_str(), msg.size());
	      csystprior_.syspriorscale[i] = 0.;
	    }

	  //Remove priors for scale variations
	  if (i >= nsyspdf)
	    {
	      /*
	      char nuispar[64];
	      sprintf (nuispar, "scale_nuisance_param_%02d", i+1 - nsyspdf);
	      cout << "Remove prior for syst " << nuispar << endl;
	      string msg = (string) "I: Remove prior for systematic " + nuispar;
	      hf_errlog_(15082401, msg.c_str(), msg.size());
	      csystprior_.syspriorscale[i] = 0.;
	      */

	      //systscal_.sysform[i] = 4; //External (Minuit Fit)
	      //systscal_.sysform[i] = 3; //Offset (Minuit Fit)
	    }
	  
	}
	  
      //Add the PDF uncertainties to the total number of systematic uncertainties
      systema_.nsys = nsysloc;
      systematicsflags_.resetcommonsyst = true;

      //double chi2tot = chi2data_theory_(2);
      string fname = outdir + "/Results_000.txt";
      fopen_(85, fname.c_str(), fname.size());
      double chi2tot = chi2data_theory_(3);
      fclose_(85);
      writefittedpoints_(); //write out fittedresults.txt file
      bool mv = system(((string)"mv " + outdir + "/fittedresults.txt " 
			+ outdir + "/fittedresults_000.txt").c_str());
	  
      char chi2c[20];
      sprintf(chi2c, "%.1f", chi2tot);
      cout << setw(15) << "chi2=" << chi2c 
	   << setw(15) << "ndf=" << cfcn_.ndfmini 
	   << endl;
    }

  //variables which stores results
  double min, deltap, deltam, chi2min;
  double central, eplus, eminus;
  vector <double> deltapi; //uncertainty decomposition
  vector <double> deltami; //uncertainty decomposition

  //decomposition
  map <int, map <int, map <double, double> > > systchi2; //map of offset chi2 values for [plus/minus][systematic uncertainties][parameters value]
  
  if (lhapdfset.size() == 0)
    {
      string msg = "S: Call to.alphasscan but lhapdfset is not specified";
      hf_errlog_(18082001, msg.c_str(), msg.size());
    }

  //Start the alphas scan
  cout << "-------------------------------------------" << endl;
  cout << "Alphas Scan" << endl;

  map <double, double> chi2; //map of parameters value and chi2 values

  //LHAPDF::initPDFSet(alphaslhapdf.c_str());
  //Dyturbo::pdfname = alphaslhapdf;
  gNode["set"] = alphaslhapdf;

  //Loop on parameters points
  int iset = 0;
  for (vector <double>::iterator vit = values.begin(); vit != values.end(); vit++)
    {
      //LHAPDF::initPDF(iset);
      //c_alphas_.alphas = LHAPDF::alphasPDF(boson_masses_.Mz);
      //clhapdf_.ilhapdfset = iset;
      //Dyturbo::pdfmember = iset;
      gNode["member"] = iset;
      evol->atConfigurationChange();
      c_alphas_.alphas = evol->getAlphaS(boson_masses_.Mz);

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
	   << setw(15) << "ndf=" << cfcn_.ndfmini 
	   << endl;
      chi2[*vit] = chi2tot;

      if (decomposition)
	decompose(systchi2, *vit);

      iset++;
    }
  fitchi2_and_store (chi2, min, deltap, deltam, chi2min, "chi2scan.txt");

  if (decomposition)
    decompose_fits(systchi2, min, deltapi, deltami);

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

  //******************************************//
  //Store PDF members for plots
  if (lhapdfprofile)
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
      getpdfunctype_as(MonteCarloPDFErr, AsymHessPDFErr, SymmHessPDFErr, evol);
      
      //Start loop on PDF sets and members
      int cset = 0; //counter on pdf variations;
      int isys = 1; //counter on PDF error index;
      for (int pdfset = 0; pdfset < sets; pdfset++)
	{
	  if (pdfset ==  0)
	    {
	      //LHAPDF::initPDFSet(lhapdfset.c_str());
	      //strcpy(clhapdf_.lhapdfset,lhapdfset.c_str());
	      gNode["set"] = lhapdfset;
	    }
	  else if (pdfset ==  1)
	    {
	      if (lhapdfvarset == "")
		continue;
	      //LHAPDF::initPDFSet(lhapdfvarset.c_str());
	      //strcpy(clhapdf_.lhapdfset,lhapdfvarset.c_str());
	      gNode["set"] = lhapdfvarset;
	    }
	  gNode["member"] = 0;
	  evol->atConfigurationChange();
	      
	  //Number of PDF members
	  //int nsets = LHAPDF::numberPDF();
	  int nsets = evol->getPropertyI("NumMembers")-1;
	  cout << "Number of PDF members for this set: " << nsets << endl;
	  totset += nsets;
	  
	  //Loop on PDF members
	  for (int iset = 0; iset <= nsets; iset++)
	    {
	      //skip central member of the VAR set
	      if (pdfset == 1 && iset == 0)
		continue;
	      
	      //Init PDF member and alphas(MZ)
	      //LHAPDF::initPDF(iset);
	      //clhapdf_.ilhapdfset = iset;
	      //c_alphas_.alphas = LHAPDF::alphasPDF(boson_masses_.Mz);
	      gNode["member"] = iset;
	      evol->atConfigurationChange();
	      c_alphas_.alphas = evol->getAlphaS(boson_masses_.Mz);
	      
	      //In VAR PDF set determine if it is a model or parametrisation variation
	      if (pdfset == 1)
		{
		  MonteCarloPDFErr = 0;
		  AsymHessPDFErr = 0;
		  SymmHessPDFErr = 0;
		  ParPDFErr = 0;
		  ModPDFErr = 0;
		  if (iset > (nsets - chi2scan_.chi2nparvar))
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
		  //fill_c_common_();
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
  cout << "Chi2 at minimum: " << chi2min << "  "  << "ndf=" << (cfcn_.ndfmini-1) << endl;
  if (decomposition)
    {
      //uncertainties decomposition
      double statp = 0;
      double statm = 0;
      double systp = 0;
      double systm = 0;
      double PDFp = 0;
      double PDFm= 0;
      double scalesp = 0;
      double scalesm = 0;
      for (int s = 0; s < systema_.nsys; s++)
	{
	  string nuislabel = string(systema_.system[s]);
	  nuislabel.erase(nuislabel.find_first_of(" "), string::npos); // trim trailing whitespaces
	  
	  //cout << nuislabel << "\t+" << deltapi[s] << "\t-" << deltami[s] << endl;

	  ofstream func;
	  func.open((outdir + "/unc_dec.txt").c_str(), ofstream::out | ofstream::app);
	  func << setprecision(6) << nuislabel << "\t" << deltapi[s] << "\t" << deltami[s] << endl;
	  func.close();

	  if (nuislabel.find("PDF_nuisance_param_") != string::npos)
	    {
	      PDFp += pow(max(max(0.,  deltapi[s]),  deltami[s]), 2);
	      PDFm += pow(max(max(0., -deltapi[s]), -deltami[s]), 2);
	    }
	  else if (nuislabel.find("scale_nuisance_param_") != string::npos)
	    {
	      scalesp += pow(max(max(0.,  deltapi[s]),  deltami[s]), 2);
	      scalesm += pow(max(max(0., -deltapi[s]), -deltami[s]), 2);
	    }
	  else if (nuislabel.find("CStat_") != string::npos || nuislabel.find("stat") != string::npos)
	    {
	      statp += pow(max(max(0.,  deltapi[s]),  deltami[s]), 2);
	      statm += pow(max(max(0., -deltapi[s]), -deltami[s]), 2);
	    }
	  else
	    {
	      systp += pow(max(max(0.,  deltapi[s]),  deltami[s]), 2);
	      systm += pow(max(max(0., -deltapi[s]), -deltami[s]), 2);
	    }
	}
      //Loop on statistical uncertainties
      int npoints = cndatapoints_.npoints;
      for (int p = 0; p < npoints; p++)
	{
	  int idx = p+systema_.nsys;
	  char statname[100];
	  sprintf(statname, "stat_%d", p);

	  string statlabel = string(statname);

	  //cout << statlabel << "\t" << deltapi[idx] << "\t" << deltami[idx] << endl;

	  ofstream func;
	  func.open((outdir + "/unc_dec.txt").c_str(), ofstream::out | ofstream::app);
	  func << setprecision(6) << statlabel << "\t" << deltapi[idx] << "\t" << deltami[idx] << endl;
	  func.close();

	  statp += pow(max(max(0.,  deltapi[idx]),  deltami[idx]), 2);
	  statm += pow(max(max(0., -deltapi[idx]), -deltami[idx]), 2);
	}      
      cout << endl;
      cout << "Uncertainty decomposition" << endl;
      cout << "Statistical" << "\t+" << sqrt(statp) << "\t-" << sqrt(statm) << endl;
      cout << "Systematic"  << "\t+" << sqrt(systp) << "\t-" << sqrt(systm) << endl;
      cout << "PDFs"        << "\t\t+" << sqrt(PDFp) << "\t-" << sqrt(PDFm) << endl;
      cout << "QCD scales"  << "\t+" << sqrt(scalesp) << "\t-" << sqrt(scalesm) << endl;
      cout << "Total (decomposed)" << "\t+" << sqrt(statp+systp+PDFp+scalesp) << "\t-" <<  sqrt(statm+systm+PDFm+scalesm) << endl;
      cout << "Total (from fit)" << "\t+" << deltap << "\t-" << deltam << endl;

      ofstream func;
      func.open((outdir + "/unc_summary.txt").c_str(), ofstream::out | ofstream::app);
      func << setprecision(6);
      func << "Uncertainty decomposition" << endl;
      func << "Statistical" << "\t+" << sqrt(statp) << "\t-" << sqrt(statm) << endl;
      func << "Systematic"  << "\t+" << sqrt(systp) << "\t-" << sqrt(systm) << endl;
      func << "PDFs"        << "\t\t+" << sqrt(PDFp) << "\t-" << sqrt(PDFm) << endl;
      func << "QCD scales"  << "\t+" << sqrt(scalesp) << "\t-" << sqrt(scalesm) << endl;
      func << "Total (decomposed)" << "\t+" << sqrt(statp+systp+PDFp+scalesp) << "\t-" <<  sqrt(statm+systm+PDFm+scalesm) << endl;
      func << "Total (from fit)" << "\t+" << deltap << "\t-" << deltam << endl;
      func.close();
    }
      
  cout << endl;
}
//#endif

void asscan::decompose(map <int, map <int, map <double, double> > > &systchi2, double value)
{
  cout << "Start uncertainty decomposition" << endl;
  int npoints = cndatapoints_.npoints;
  /********************** Technique 1 *******************************
   //remove one-by-one each uncertainty from the chi2 and recalculate chi2 (should apply shift)
	  
   //loop on systematic uncertainties
	  int totsyst = systema_.nsys;
	  for (int s = 0; s < totsyst; s++)
	    {
	      double savebetaasym[npoints][2][totsyst];
	      double savebeta[npoints][totsyst];
	      double saveomega[npoints][totsyst];
	      double savetheo[npoints];
	      double savedata[npoints];
	      for (int p = 0; p < npoints; p++)
		{
		//save current uncertainty
		  savebetaasym[p][0][s] = systasym_.betaasym[p][0][s];
		  savebetaasym[p][1][s] = systasym_.betaasym[p][1][s];
		  savebeta[p][s] = systema_.beta[p][s];
		  saveomega[p][s] = systasym_.omega[p][s];

		  //save current theory
		  savetheo[p] = c_theo_.theo[p];
		  //save current data
		  savedata[p] = indata2_.daten_[p];

		  //remove uncertainty
		  systasym_.betaasym[p][0][s] = 0;
		  systasym_.betaasym[p][1][s] = 0;
		  systema_.beta[p][s] = 0;
		  systasym_.omega[p][s] = 0;

		  //shift the data (should know the shift at the chi2 minimum, need two iterations for that)
		  //indata2_.daten_[p] = indata2_.daten_[p]
		  //+ systexport_.sysshift_[s]*systexport_.scgamma_[p][s]
		  //+ systexport_.sysshift_[s]*systexport_.sysshift_[s]*systexport_.s.omega[p][s];

		  //offset the data
		  //indata2_.daten_[p] = indata2_.daten_[p] + systexport_.scgamma_[p][s] + systexport_.s.omega[p][s];
		  //indata2_.daten_[p] = indata2_.daten_[p] + savebeta[p][s]*savetheo[p] + saveomega[p][s]*savetheo[p];
		} //end loop on points

		//calculate chi2
	      systchi2[s][value] = chi2data_theory_(2);

	      //char chi2c[20];
	      //sprintf(chi2c, "%.8f", systchi2[s][value]);
	      //cout << setw(15) << (label + "=") << value
	      //	   << setw(6) << "syst" << setw(4) << s // << setw(20) << systema_.system[s]
	      //	   << setw(15) << "chi2=" << chi2c 
	      //	   << setw(15) << "ndf=" << cfcn_.ndfmini 
	      //     << endl;

	      for (int p = 0; p < npoints; p++)
		{
		//restore uncertainty
		  systasym_.betaasym[p][0][s] = savebetaasym[p][0][s];
		  systasym_.betaasym[p][1][s] = savebetaasym[p][1][s];
		  systema_.beta[p][s] = savebeta[p][s];
		  systasym_.omega[p][s] = saveomega[p][s];


		  ////restore theory
		  //c_theo_.theo[p] = savetheo[p];
		  ////restore data
		  //indata2_.daten_[p] = savedata[p];

		} //end loop on points
		  
	    }//end loop on uncertainties
  */
	  
  /********************** Technique 2 *******************************/
  //offset the data one by one

  int totsyst = systema_.nsys;
  //theory and data
  double savetheo[npoints];
  double savedata[npoints];
  //systematic uncertainties
  double savebetaasym[npoints][2][totsyst];
  double savebeta[npoints][totsyst];
  double saveomega[npoints][totsyst];
  double savescgamma[npoints][totsyst];
  double savescomega[npoints][totsyst];
  //statistical uncertainties
  double savestatpoi[npoints];
  double savestatconst[npoints];
  double saveuncorpoi[npoints];
  double saveuncorconst[npoints];
	      
  for (int p = 0; p < npoints; p++)
    {
      for (int s = 0; s < systema_.nsys; s++)
	{
	  //save current uncertainty
	  savebetaasym[p][0][s] = systasym_.betaasym[p][0][s];
	  savebetaasym[p][1][s] = systasym_.betaasym[p][1][s];
	  savebeta[p][s] = systema_.beta[p][s];
	  saveomega[p][s] = systasym_.omega[p][s];

	  //save current scaled gamma and omega
	  savescgamma[p][s] = systexport_.scgamma_[p][s];
	  savescomega[p][s] = systexport_.scomega_[p][s];
	}
	      
      //save current theory
      savetheo[p] = c_theo_.theo[p];
      //save current data
      savedata[p] = indata2_.daten_[p];
	      
      //save current stat uncertainties
      savestatpoi[p] = cuncerrors_.e_stat_poisson[p];
      savestatconst[p] = cuncerrors_.e_stat_const[p];
      saveuncorpoi[p] = cuncerrors_.e_uncor_poisson[p];
      saveuncorconst[p] = cuncerrors_.e_uncor_const[p];
    }
	  
  //loop on systematic uncertainties
  for (int s = 0; s < totsyst; s++)
    for (int sign = 0; sign < 2; sign++)
      {
	//offset the data
	for (int p = 0; p < npoints; p++)
	  {
	    //offset the data by one standard deviation of the current systematic uncertainty
	    if (sign == 0)
	      {
		indata2_.daten_[p] = savedata[p] + savescgamma[p][s] + savescomega[p][s];                     //These are the Gamma eventually scaled to the theory
		//if (systscal_.sysscalingtype[s] == 0) //data-like uncertainty, scaled to the data
		//  indata2_.daten_[p] = savedata[p] + savebeta[p][s]*savedata[p] + saveomega[p][s]*savedata[p];    //These are the original unscaled Gamma
		//else if (systscal_.sysscalingtype[s] == 1) //theory-like uncertainty, scaled to the theory
		//  indata2_.daten_[p] = savedata[p] + savebeta[p][s]*savetheo[p] + saveomega[p][s]*savetheo[p];  //These are the original Gamma scaled to the theory
		//else if (systscal_.sysscalingtype[s] == 2) //poissonian
		//  indata2_.daten_[p] = savedata[p] + savebeta[p][s]*sqrt(savetheo[p]*savedata[p]) + saveomega[p][s]*sqrt(savetheo[p]*savedata[p]);
	      }
	    else
	      {
		indata2_.daten_[p] = savedata[p] - savescgamma[p][s] + savescomega[p][s];                     //These are the Gamma eventually scaled to the theory
		//if (systscal_.sysscalingtype[s] == 0) //data-like uncertainty, scaled to the data
		//  indata2_.daten_[p] = savedata[p] - savebeta[p][s]*savedata[p] + saveomega[p][s]*savedata[p];    //These are the original unscaled Gamma
		//else if (systscal_.sysscalingtype[s] == 1) //theory-like uncertainty, scaled to the theory
		//  indata2_.daten_[p] = savedata[p] - savebeta[p][s]*savetheo[p] + saveomega[p][s]*savetheo[p];  //These are the original Gamma scaled to the theory
		//else if (systscal_.sysscalingtype[s] == 2) //poissonian
		//  indata2_.daten_[p] = savedata[p] - savebeta[p][s]*sqrt(savetheo[p]*savedata[p]) + saveomega[p][s]*sqrt(savetheo[p]*savedata[p]);
	      }
	    //cout << p << " daten " << indata2_.daten_[p] << " savedata " <<  savedata[p] << " savegamma " <<  savescgamma[p][s] << " saveomega " << savescomega[p][s] << endl;
	  }

	//recompute beta and omega, and uncorrelated uncertainties, so that absolute errors are unchanged
	for (int p = 0; p < npoints; p++)
	  {
	    if (indata2_.daten_[p] == 0) continue;

	    for (int s2 = 0; s2 < totsyst; s2++)
	      {
		double fac = 1;
		if (systscal_.sysscalingtype[s2] == 0) //additive uncertainty
		  fac = savedata[p]/indata2_.daten_[p]; 
		else if (systscal_.sysscalingtype[s2] == 1) //multiplicative uncertainty
		  fac = 1.;
		else if (systscal_.sysscalingtype[s2] == 2) //poissonian uncertainty
		  fac = sqrt(savetheo[p]*savedata[p])/sqrt(indata2_.daten_[p]*savetheo[p]);
			
		systasym_.betaasym[p][0][s2] = savebetaasym[p][0][s2]*fac;
		systasym_.betaasym[p][1][s2] = savebetaasym[p][1][s2]*fac;
		systema_.beta[p][s2] = savebeta[p][s2]*fac;
		systasym_.omega[p][s2] = saveomega[p][s2]*fac;
	      }
	    cuncerrors_.e_stat_const[p]    = savestatconst[p]*savedata[p]/indata2_.daten_[p];
	    cuncerrors_.e_uncor_const[p]   = saveuncorconst[p]*savedata[p]/indata2_.daten_[p];
	    if (cuncerrors_.e_stat_poisson[p] > 0 && indata2_.daten_[p]*savetheo[p] > 0)
	      cuncerrors_.e_stat_poisson[p]  = savestatpoi[p]*sqrt(savetheo[p]*savedata[p])/sqrt(indata2_.daten_[p]*savetheo[p]);
	    if (cuncerrors_.e_uncor_poisson[p] > 0 && indata2_.daten_[p]*savetheo[p] > 0)
	      cuncerrors_.e_uncor_poisson[p] = saveuncorpoi[p]*sqrt(savetheo[p]*savedata[p])/sqrt(indata2_.daten_[p]*savetheo[p]);
	  } //end loop on points
		
	//calculate chi2
	//systchi2[sign][s][value] = chi2data_theory_(2);

	int iflag = 2;
	int n0 = npoints;
	double fchi2 = 0.;
	double rsys[totsyst];
	double ersys[totsyst];
	for (int i = 0; i < totsyst; i++)
	  {
	    rsys[i] = 0.;
	    ersys[i] = 0.;
	  }
	double pchi2[NSET_C];
	double fcorchi2 = 0.;
	getnewchisquare_(iflag,n0,fchi2,rsys,ersys,pchi2,fcorchi2);
	systchi2[sign][s][value] = fchi2;
		  
	//char chi2c[20];
	//sprintf(chi2c, "%.8f", systchi2[sign][s][value]);
	//cout << setw(15) << value
	//     << setw(6) << "syst" << setw(4) << s // << setw(20) << systema_.system[s]
	//     << setw(15) << "chi2=" << chi2c 
	//     << setw(15) << "ndf=" << cfcn_.ndfmini
	//     << endl;
		
	for (int p = 0; p < npoints; p++)
	  {
	    //restore uncertainty
	    for (int s2 = 0; s2 < totsyst; s2++)
	      {
		systasym_.betaasym[p][0][s2] = savebetaasym[p][0][s2];
		systasym_.betaasym[p][1][s2] = savebetaasym[p][1][s2];
		systema_.beta[p][s2] = savebeta[p][s2];
		systasym_.omega[p][s2] = saveomega[p][s2];
	      }

	    //restore theory
	    c_theo_.theo[p] = savetheo[p];
	    //restore data
	    indata2_.daten_[p] = savedata[p];
		    
	    //restore stat uncertainties
	    cuncerrors_.e_stat_poisson[p] = savestatpoi[p];
	    cuncerrors_.e_stat_const[p] = savestatconst[p];
	    cuncerrors_.e_uncor_poisson[p] = saveuncorpoi[p];
	    cuncerrors_.e_uncor_const[p] = saveuncorconst[p];
	  } //end loop on points
      }//end loop on systematic uncertainties and plus/minus
	
  //loop on statistical uncertainties (also loop on points)
  for (int p = 0; p < npoints; p++)
    for (int sign = 0; sign < 2; sign++)
      {
	//offset the data
	if (sign == 0)
	  indata2_.daten_[p] = savedata[p] + systexport_.scerrors_[p];
	else
	  indata2_.daten_[p] = savedata[p] - systexport_.scerrors_[p];

	//recompute beta and omega, and stat uncertainties, so that scaled errors are invariant
	//if (indata2_.daten_[p] == 0) continue;

	for (int s2 = 0; s2 < totsyst; s2++)
	  {
	    double fac = 1.;
	    if (systscal_.sysscalingtype[s2] == 0) //additive uncertainty
	      fac = savedata[p]/indata2_.daten_[p];
	    else if (systscal_.sysscalingtype[s2] == 1)
	      fac = 1.;
	    else if (systscal_.sysscalingtype[s2] == 2)
	      fac = sqrt(savetheo[p]*savedata[p])/sqrt(indata2_.daten_[p]*savetheo[p]);
			
	    systasym_.betaasym[p][0][s2] = savebetaasym[p][0][s2]*fac;
	    systasym_.betaasym[p][1][s2] = savebetaasym[p][1][s2]*fac;
	    systema_.beta[p][s2] = savebeta[p][s2]*fac;
	    systasym_.omega[p][s2] = saveomega[p][s2]*fac;
	  }
	cuncerrors_.e_stat_const[p]    = savestatconst[p]*savedata[p]/indata2_.daten_[p];
	cuncerrors_.e_uncor_const[p]   = saveuncorconst[p]*savedata[p]/indata2_.daten_[p];
	if (cuncerrors_.e_stat_poisson[p] > 0 && indata2_.daten_[p]*savetheo[p] > 0)
	  cuncerrors_.e_stat_poisson[p]  = savestatpoi[p]*sqrt(savetheo[p]*savedata[p])/sqrt(indata2_.daten_[p]*savetheo[p]);
	if (cuncerrors_.e_uncor_poisson[p] > 0 && indata2_.daten_[p]*savetheo[p] > 0)
	  cuncerrors_.e_uncor_poisson[p] = saveuncorpoi[p]*sqrt(savetheo[p]*savedata[p])/sqrt(indata2_.daten_[p]*savetheo[p]);

	//calculate chi2
	//systchi2[sign][p+totsyst][value] = chi2data_theory_(2);

	int iflag = 2;
	int n0 = npoints;
	double fchi2 = 0.;
	double rsys[totsyst];
	double ersys[totsyst];
	for (int i = 0; i < totsyst; i++)
	  {
	    rsys[i] = 0.;
	    ersys[i] = 0.;
	  }
	double pchi2[NSET_C];
	double fcorchi2 = 0.;
	getnewchisquare_(iflag,n0,fchi2,rsys,ersys,pchi2,fcorchi2);
	systchi2[sign][p+totsyst][value] = fchi2;
	
	//char chi2c[20];
	//sprintf(chi2c, "%.8f", systchi2[sign][p+totsyst][value]);
	//cout << setw(15) << value
	//     << setw(6) << "stat" << setw(4) << p // << setw(20) << systema_.system[s]
	//     << setw(15) << "chi2=" << chi2c 
	//     << setw(15) << "ndf=" << cfcn_.ndfmini 
	//     << endl;
		  
	//restore uncertainty
	for (int s = 0; s < totsyst; s++)
	  {
	    systasym_.betaasym[p][0][s] = savebetaasym[p][0][s];
	    systasym_.betaasym[p][1][s] = savebetaasym[p][1][s];
	    systema_.beta[p][s] = savebeta[p][s];
	    systasym_.omega[p][s] = saveomega[p][s];
	  }

	//restore theory
	c_theo_.theo[p] = savetheo[p];
	//restore data
	indata2_.daten_[p] = savedata[p];

	//restore stat uncertainties
	cuncerrors_.e_stat_poisson[p] = savestatpoi[p];
	cuncerrors_.e_stat_const[p] = savestatconst[p];
	cuncerrors_.e_uncor_poisson[p] = saveuncorpoi[p];
	cuncerrors_.e_uncor_const[p] = saveuncorconst[p];
      } //end loop on stat uncertainties and plus/minus
  /***************************************************/
}

void asscan::decompose_fits(map <int, map <int, map <double, double> > > systchi2, double min, vector <double> &deltapi, vector <double> &deltami)
{
  int npoints = cndatapoints_.npoints;
  /***************** Technique 1 ******************/
  /*
  //Loop on uncertainties
  double deltap2_tot = 0;
  double deltam2_tot = 0;
  for (int s = 0; s < systema_.nsys; s++)
  {
  double min_i, deltap_i, deltam_i, chi2min_i;
  char chi2name[100];
  sprintf(chi2name, "chi2scan_syst_%d.txt", s);
  fitchi2_and_store (systchi2[s], min_i, deltap_i, deltam_i, chi2min_i, chi2name);
  deltap2_tot += max(0.,deltap*deltap-deltap_i*deltap_i);
  deltam2_tot += max(0.,deltam*deltam-deltam_i*deltam_i);

  deltapi.push_back(sqrt(max(0.,deltap*deltap-deltap_i*deltap_i)));
  deltami.push_back(sqrt(max(0.,deltam*deltam-deltam_i*deltam_i)));
  }
  */

  /***************** Technique 2 ******************/
  //Loop on systematic uncertainties
  for (int s = 0; s < systema_.nsys; s++)
    {
      double min_i_p, min_i_m, deltap_i, deltam_i, chi2min_i;
      char chi2name[100];
      sprintf(chi2name, "chi2scan_syst_%d_p.txt", s);
      fitchi2_and_store (systchi2[0][s], min_i_p, deltap_i, deltam_i, chi2min_i, chi2name);

      sprintf(chi2name, "chi2scan_syst_%d_m.txt", s);
      fitchi2_and_store (systchi2[1][s], min_i_m, deltap_i, deltam_i, chi2min_i, chi2name);

      deltapi.push_back(min_i_p - min);
      deltami.push_back(min_i_m - min);

      //deltapi.push_back(max(max(0., min_i_p - min), min_i_m - min));
      //deltami.push_back(max(max(0., min - min_i_p), min - min_i_m));
    }
  //Loop on statistical uncertainties
  for (int p = 0; p < npoints; p++)
    {
      int idx = p+systema_.nsys;
      double min_i_p, min_i_m, deltap_i, deltam_i, chi2min_i;
      char chi2name[100];
      sprintf(chi2name, "chi2scan_syst_%d_p.txt", idx);
      fitchi2_and_store (systchi2[0][idx], min_i_p, deltap_i, deltam_i, chi2min_i, chi2name);

      sprintf(chi2name, "chi2scan_syst_%d_m.txt", idx);
      fitchi2_and_store (systchi2[1][idx], min_i_m, deltap_i, deltam_i, chi2min_i, chi2name);

      deltapi.push_back(min_i_p - min);
      deltami.push_back(min_i_m - min);

      //deltapi.push_back(max(max(0., min_i_p - min), min_i_m - min));
      //deltami.push_back(max(max(0., min - min_i_p), min - min_i_m));
    }
  /***********************************/
}
