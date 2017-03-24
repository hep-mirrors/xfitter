#include "xfitter_cpp.h"

#include "pdferrors.h"
#include "dimensions.h"

#ifdef LHAPDF_ENABLED
#include <LHAPDF/LHAPDF.h>
#endif
#define LHAPDF6 6

#include "TheorEval.h"

#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <math.h>
#include <numeric>
#include <algorithm>

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


vector <double> mcweights(vector<double> const& chi2, int ndata, bool GK_method)
{
  vector<double> w;
  const int nrep = chi2.size();
  
  vector<double> logw;  // Calculate Ln(wnn) (nn - not normalised) -> use repnums as index for unweighted PDFs
  for (vector <double>::const_iterator c = chi2.begin(); c != chi2.end(); c++) {
    if ( (GK_method) == false) {
      logw.push_back( - (*c)/(2.0) +( (((double) ndata)-1.0)/2.)*log(*c)); //Bayesian
    }
    if ( (GK_method) == true)  logw.push_back(- (*c)/(2.0)); //Giele-Keller
  }

  // Get maximum value of log(w)
  double exp_avg = *(max_element(logw.begin(), logw.end()));
  //exp_avg =0;//only needed to normalise

  // Calculate weights
  for (size_t i=0 ;i<nrep; i++) {
    w.push_back(exp(logw[i] - exp_avg ));
  }

  //Drop any weights smaller than 1e-12
  double wtot = accumulate(w.begin(),w.end(),0.0); 

  for (size_t i=0; i<w.size(); i++) {
    if ((w[i]*(nrep/wtot)) < 1e-12)
      w[i]=0;
  }

     wtot=accumulate(w.begin(),w.end(),0.0);  // Normalise weights so Sum(weights)=N

     for (size_t i=0;i<w.size();i++){
     	w[i]*=(nrep/wtot); 
    }

   return w;
}


//return error if LHAPDF is not enabled
#ifndef LHAPDF_ENABLED
void get_lhapdferrors_()
{
  string msg = "S: Call to LHAPDFErrors but LHAPDF is not enabled. Run ./configure --enable-lhapdf and link the executable";
  hf_errlog_(14060101, msg.c_str(), msg.size());
}

void set_verbosity_(int& i)
{
  string msg = "S: Call to LHAPDF set_verbosity but LHAPDF is not enabled. Run ./configure --enable-lhapdf and link the executable";
  hf_errlog_(17030201, msg.c_str(), msg.size());
}

#else

void set_verbosity_(int& i)
{
  LHAPDF::setVerbosity(i);
}

    

void get_lhapdferrors_()
{
  cout << endl << endl << endl;
  cout << "  -----------------------" << endl;
  cout << "  Start LHAPDFErrors chi2 evaluation" << endl;
  cout << "  -----------------------" << endl;
  cout << endl << endl << endl;

  
  //Start program
  int npoints = cndatapoints_.npoints_;
  int nsysloc = systema_.nsys_;

  //Number of systematic uncertainties
  //cout << " NSYST = " << nsysloc << "\n";

  int MonteCarloPDFErr = 0;
  int AsymHessPDFErr = 0;
  int SymmHessPDFErr = 0;
  int ModPDFErr = 0;
  int ParPDFErr = 0;

  map <int, point> pointsmap;
  vector <double> chi2;

  string lhapdfset = string(clhapdf_.lhapdfset_, 128);
  lhapdfset = lhapdfset.erase(lhapdfset.find_last_not_of(" ")+1, string::npos);
  string lhapdfvarset = string(clhapdf_.lhapdfvarset_, 128);
  lhapdfvarset = lhapdfvarset.erase(lhapdfvarset.find_last_not_of(" ")+1, string::npos);
  bool pdfprofile = clhapdf_.lhapdfprofile_;
  bool scaleprofile = clhapdf_.lhascaleprofile_;

  string outdirname = string(coutdirname_.outdirname_, 256);
  outdirname = outdirname.erase(outdirname.find_last_not_of(" ")+1, string::npos);

  /*****************************************************/
  // Central PDF
  cout << "-------------------------------------------" << endl;
  cout << "Chi2 test on central prediction:" << endl;
  
  LHAPDF::initPDFSet(lhapdfset.c_str());

  int central_pdfmember;
  if (pdfprofile)
    central_pdfmember = 0;
  else
    central_pdfmember = clhapdf_.ilhapdfset_;
  clhapdf_.ilhapdfset_ = central_pdfmember;
  LHAPDF::initPDF(clhapdf_.ilhapdfset_);

  //set alphas from LHAPDF
  c_alphas_.alphas_ = LHAPDF::alphasPDF(boson_masses_.mz_);

#if LHAPDF_FAMILY == LHAPDF6
  //set also mc, mb and mt from LHAPDF
  steering_.hf_mass_[0] = LHAPDF::getThreshold(4);
  steering_.hf_mass_[1] = LHAPDF::getThreshold(5);
  steering_.hf_mass_[2] = LHAPDF::getThreshold(6);
#endif

  //Evaluate chi2 for the central PDF member
  double chi2tot;
  //  mc_method_();
  chi2tot = chi2data_theory_(1); //Needed for initialiion
  string fname = outdirname + "/Results.txt";
  fopen_(85, fname.c_str(), fname.size());
  chi2tot = chi2data_theory_(3);
  fclose_(85);
  bool cp = system(((string)"cp " + outdirname + "/Results.txt " + outdirname + "/Results_00.txt").c_str());

  char chi2c[500];
  sprintf(chi2c, "%.2f", chi2tot);

  cout << setw(20) << "PDF set number: " << setw(5) << central_pdfmember
       << setw(15) << "chi2/ndf = " << chi2c << "/" << cfcn_.ndfmini_
       << endl;

  //Write results of chi2 calculation in the fittedresults.txt file
  char tag[10];
  sprintf (tag, "_%04d", central_pdfmember);
  writefittedpoints_(); //write out 
  cp = system(((string)"cp " + outdirname + "/fittedresults.txt " 
		    + outdirname + "/fittedresults.txt_set" + tag).c_str());
  //  writefittedpoints_();
      
  //Store PDF member
  string filename = outdirname + "/pdfs_q2val_";
  filename += " ";
  store_pdfs_(filename.c_str(), filename.size());
      
  //Save the input PDFs in LHAPDF6 format
  if ( strstr(coutdirname_.lhapdf6outdir_ ,"xfitter_pdf") != NULL) { 
    // Overwrite default name with the PDF name
    string msg = (string) "I: Overwrite the lhapdf6 output dir name  with input PDF name, "+  clhapdf_.lhapdfset_;
    hf_errlog_(15051301,msg.c_str(), msg.size());
    strncpy(coutdirname_.lhapdf6outdir_, clhapdf_.lhapdfset_, 256);
  }  
  fill_c_common_(); //set the PDF name and the HF masses in the common block, so they are written out correctly in the LHAPDF6 output
  print_lhapdf6_();
  //save_data_lhapdf6_(pdfmember);

  //Store the central theory prediction
  for (int i = 0; i < npoints; i++)
    pointsmap[i].thc = c_theo_.theo_[i];
  /*****************************************************/

#if LHAPDF_FAMILY == LHAPDF6
  //Reduce LHAPDF verbosity
  LHAPDF::Info& cfg = LHAPDF::getConfig();
  cfg.set_entry("Verbosity", 0);
#endif
  
  
  /*****************************************************/
  // PDF variations
  int cset = 1; //global counter on total pdf and scale variations;
  if (pdfprofile)
    {
      cout << "-------------------------------------------" << endl;
      cout << "Chi2 test on PDF variations" << endl;
      getpdfunctype_heraf_(MonteCarloPDFErr, AsymHessPDFErr, SymmHessPDFErr, lhapdfset.c_str(), lhapdfset.size());
      string msg = "";
      if (MonteCarloPDFErr)
	msg = (string) "I: Use Monte Carlo errors approach for: " + lhapdfset;
      else if (SymmHessPDFErr)
	msg = (string) "I: Use symmetric hessian errors approach for: " + lhapdfset;
      else if(AsymHessPDFErr)
	msg = (string) "I: Use asymmetric hessian errors approach for: " + lhapdfset;
      else
	msg = (string) "S: Error determining PDF errors approach for: " + lhapdfset;
      hf_errlog_(14012701, msg.c_str(), msg.size());

      if (clhapdf_.scale68_)
	{
	  msg = (string) "I: Scale Errors from 90cl to 68cl for: " + lhapdfset;
	  hf_errlog_(23051401, msg.c_str(), msg.size());
	}

      MonteCarloPDFErr = 0;
      AsymHessPDFErr = 0;
      SymmHessPDFErr = 0;
      ParPDFErr = 0;
      ModPDFErr = 0;

      //Loop on PDF members, and on EIG and VAR pdf sets if available
      int isys = 1; //counter on PDF error index;
      for (int pdfset = 0; pdfset < 2; pdfset++)
	{
	  if (pdfset ==  0)
	    {
	      LHAPDF::initPDFSet(lhapdfset.c_str());
	      //	  LHAPDF::initPDFSetByName(lhapdfset.c_str());
	      cout << "PDF set: " << lhapdfset;
	      getpdfunctype_heraf_(MonteCarloPDFErr, AsymHessPDFErr, SymmHessPDFErr, lhapdfset.c_str(), lhapdfset.size());
	    }
	  else if (pdfset ==  1)
	    {
	      if (lhapdfvarset == "")
		continue;
	      LHAPDF::initPDFSet(lhapdfvarset.c_str());
	    }
	  //Number of PDF members
	  int nsets = LHAPDF::numberPDF();
      
	  cout << ", number of PDF variations for this set: " << nsets << endl;
	  if (nsets == 1) //Needed for PDF sets with only the central member
	    nsets = 0;

	  //VAR set info message
	  if (pdfset ==  1)
	    {
	      char mod[10];
	      char par[10];
	      sprintf (par, "%d", clhapdf_.nparvar_);
	      sprintf (mod, "%d", nsets - clhapdf_.nparvar_);
	      msg = (string) "I: Add " + mod + " model and " + par + " parametrisation uncertainties from " + lhapdfvarset;
	      hf_errlog_(14012702, msg.c_str(), msg.size());
	    }

	  for (int iset = 1; iset <= nsets; iset++)
	    {
	      //skip central member of the VAR set
	      if (pdfset == 1 && iset == 0)
		continue;
	      
	      //Init PDF member and alphas(MZ)
	      LHAPDF::initPDF(iset);
	      clhapdf_.ilhapdfset_ = iset;
	      c_alphas_.alphas_ = LHAPDF::alphasPDF(boson_masses_.mz_);
	  
#if LHAPDF_FAMILY == LHAPDF6
	      //set also mc, mb and mt from LHAPDF
	      steering_.hf_mass_[0] = LHAPDF::getThreshold(4);
	      steering_.hf_mass_[1] = LHAPDF::getThreshold(5);
	      steering_.hf_mass_[2] = LHAPDF::getThreshold(6);
#endif

	      //set the HF masses in the common block, so they are written out correctly in the LHAPDF6 output
	      fill_c_common_();

	      //In VAR PDF set determine if it is a model or parametrisation variation
	      if (pdfset == 1)
		{
		  MonteCarloPDFErr = 0;
		  AsymHessPDFErr = 0;
		  SymmHessPDFErr = 0;
		  ParPDFErr = 0;
		  ModPDFErr = 0;
		  if (iset > (nsets - clhapdf_.nparvar_))
		    ParPDFErr = true;
		  else
		    ModPDFErr = true;
		}

	      //convert the PDF member index to a 4 characters string
	      char tag[10];
	      sprintf (tag, "%02d", cset);

	      //Evaluate chi2 for the current PDF member
	      double chi2tot;
	      /*
	      string fname = outdirname + "/Results_" + tag + ".txt";
	      fopen_(85, fname.c_str(), fname.size());
	      chi2tot = chi2data_theory_(3);
	      fclose_(85);
	      */
	      chi2tot = chi2data_theory_(2);

	      char chi2c[500];
	      sprintf(chi2c, "%.2f", chi2tot);

	      cout << setw(20) << "PDF set number: " << setw(5) << iset 
		   << setw(15) << "chi2/ndf = " << chi2c << "/" << cfcn_.ndfmini_
		   << endl;

	      chi2.push_back(chi2tot) ;
  
	      //Write results of chi2 calculation
	      sprintf (tag, "%04d", cset);
	      writefittedpoints_(); //write out fittedresults.txt file
	      bool mv = system(((string)"mv " + outdirname + "/fittedresults.txt " 
				+ outdirname + "/fittedresults.txt_set_" + tag).c_str());

	      //Store PDF members
	      string filename = outdirname + "/pdfs_q2val_";
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
	      
	      filename += " ";

	      store_pdfs_(filename.c_str(), filename.size());
	      
	      //Save the input PDFs in LHAPDF6 format
	      save_data_lhapdf6_(iset);

	      //Store all PDF variations
	      for (int i = 0; i < npoints; i++)
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
	  
	      if (MonteCarloPDFErr || SymmHessPDFErr || ParPDFErr)
		isys++;
	      else if ((iset%2) == 0) //set the same index for Up and Down variation of asymmetric PDF errors
		isys++;
	      cset++;
	    } //End of loop on PDF members
	} //End loop on eig var PDF sets

      // Reset the central PDF set:
      LHAPDF::initPDFSet(lhapdfset.c_str());
      LHAPDF::initPDF(central_pdfmember);
      clhapdf_.ilhapdfset_ = central_pdfmember;
      c_alphas_.alphas_ = LHAPDF::alphasPDF(boson_masses_.mz_);
#if LHAPDF_FAMILY == LHAPDF6
      steering_.hf_mass_[0] = LHAPDF::getThreshold(4);
      steering_.hf_mass_[1] = LHAPDF::getThreshold(5);
      steering_.hf_mass_[2] = LHAPDF::getThreshold(6);
#endif

    }//End of pdfprofile

  /*************************************************/
  //In QCD scales profiling mode, add two nuisance parameters for the renormalisation and factorisation scales
  if (scaleprofile)
    {
      cout << "-------------------------------------------" << endl;
      cout << "Chi2 test on QCD scale variations" << endl;
      
      //Set the central PDF
      LHAPDF::initPDFSet(lhapdfset.c_str());
      LHAPDF::initPDF(central_pdfmember);
      clhapdf_.ilhapdfset_ = central_pdfmember;
      c_alphas_.alphas_ = LHAPDF::alphasPDF(boson_masses_.mz_);
#if LHAPDF_FAMILY == LHAPDF6
      steering_.hf_mass_[0] = LHAPDF::getThreshold(4);
      steering_.hf_mass_[1] = LHAPDF::getThreshold(5);
      steering_.hf_mass_[2] = LHAPDF::getThreshold(6);
#endif

      //Store current values of scales
      map <int, int> iordmap;
      map <int, double> murmap;
      map <int, double> mufmap;

      for (map <int, TheorEval* >::iterator tit = gTEmap.begin(); tit != gTEmap.end(); tit++)
	{
	  int iord;
	  double mur0;
	  double muf0;
	  double mures0;
	  tit->second->GetOrdScales(iord, mur0, muf0);
	  iordmap[tit->first] = iord;
	  murmap[tit->first] = mur0;
	  mufmap[tit->first] = muf0;
	}

      //Apply a factor for scale variations
      double factor = 2.;
      double chi2tot;

      //mur*2
      for (map <int, TheorEval* >::iterator tit = gTEmap.begin(); tit != gTEmap.end(); tit++)
	tit->second->SetOrdScales(iordmap[tit->first], murmap[tit->first]*factor, mufmap[tit->first]);
      chi2tot = chi2data_theory_(2);
      char chi2c[500];
      sprintf(chi2c, "%.2f", chi2tot);
      cout << setw(20) << "mur = 2.0: "
	   << setw(15) << "chi2/ndf = " << chi2c << "/" << cfcn_.ndfmini_
	   << endl;
      char tag[10]; sprintf (tag, "_%04d", cset);
      writefittedpoints_(); //write out fittedresults.txt file
      bool mv = system(((string)"mv " + outdirname + "/fittedresults.txt " + outdirname + "/fittedresults.txt_scale" + tag).c_str());
      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
	pointsmap[i].th_scale_p.push_back(c_theo_.theo_[i]);
      cset++;

      //mur*0.5
      for (map <int, TheorEval* >::iterator tit = gTEmap.begin(); tit != gTEmap.end(); tit++)
	tit->second->SetOrdScales(iordmap[tit->first], murmap[tit->first]/factor, mufmap[tit->first]);
      chi2tot = chi2data_theory_(2);
      chi2c[500];
      sprintf(chi2c, "%.2f", chi2tot);
      cout << setw(20) << "mur = 0.5: "
	   << setw(15) << "chi2/ndf = " << chi2c << "/" << cfcn_.ndfmini_
	   << endl;
      tag[10]; sprintf (tag, "_%04d", cset);
      writefittedpoints_(); //write out fittedresults.txt file
      mv = system(((string)"mv " + outdirname + "/fittedresults.txt " + outdirname + "/fittedresults.txt_scale" + tag).c_str());
      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
	pointsmap[i].th_scale_m.push_back(c_theo_.theo_[i]);
      cset++;

      //muf*2
      for (map <int, TheorEval* >::iterator tit = gTEmap.begin(); tit != gTEmap.end(); tit++)
	tit->second->SetOrdScales(iordmap[tit->first], murmap[tit->first], mufmap[tit->first]*factor);
      chi2tot = chi2data_theory_(2);
      chi2c[500];
      sprintf(chi2c, "%.2f", chi2tot);
      cout << setw(20) << "muf = 2.0: "
	   << setw(15) << "chi2/ndf = " << chi2c << "/" << cfcn_.ndfmini_
	   << endl;
      tag[10]; sprintf (tag, "_%04d", cset);
      writefittedpoints_(); //write out fittedresults.txt file
      mv = system(((string)"mv " + outdirname + "/fittedresults.txt " + outdirname + "/fittedresults.txt_scale" + tag).c_str());
      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
	pointsmap[i].th_scale_p.push_back(c_theo_.theo_[i]);
      cset++;

      //muf*0.5
      for (map <int, TheorEval* >::iterator tit = gTEmap.begin(); tit != gTEmap.end(); tit++)
	tit->second->SetOrdScales(iordmap[tit->first], murmap[tit->first], mufmap[tit->first]/factor);
      chi2tot = chi2data_theory_(2);
      chi2c[500];
      sprintf(chi2c, "%.2f", chi2tot);
      cout << setw(20) << "muf = 0.5: "
	   << setw(15) << "chi2/ndf = " << chi2c << "/" << cfcn_.ndfmini_
	   << endl;
      tag[10]; sprintf (tag, "_%04d", cset);
      writefittedpoints_(); //write out fittedresults.txt file
      mv = system(((string)"mv " + outdirname + "/fittedresults.txt " + outdirname + "/fittedresults.txt_scale" + tag).c_str());
      for (int i = 0; i < npoints; i++) //Store the scale variation for each data point
	pointsmap[i].th_scale_m.push_back(c_theo_.theo_[i]);
      cset++;

      //restore nominal scale
      for (map <int, TheorEval* >::iterator tit = gTEmap.begin(); tit != gTEmap.end(); tit++)
	tit->second->SetOrdScales(iordmap[tit->first], murmap[tit->first], mufmap[tit->first]);
    }
  /*************************************************/


  //Monte Carlo replica PDF uncertainties
  int totmc = pointsmap.begin()->second.th_mc.size();;
  if (totmc > 0)
    {
      //Store weights for reweighting
      vector <double> weights = mcweights(chi2, npoints, false);
 
      ofstream chi2wf((outdirname + "/pdf_BAYweights.dat").c_str());

      vector <double>::iterator ic = chi2.begin();
      vector <double>::iterator iw = weights.begin();

      string lhapdfsetname=lhapdfset;
      if (lhapdfsetname.find(".LHgrid")!=string::npos)   lhapdfsetname.erase(lhapdfsetname.find(".LHgrid"));

      chi2wf << "LHAPDF set=   " << lhapdfsetname<<endl;
      chi2wf << "Reweight method=   BAYESIAN"<<endl;
      chi2wf << "ndata=   " << npoints <<endl;

      chi2wf << totmc << endl;

      for (; ic != chi2.end(); ic++, iw++)
	chi2wf << (ic - chi2.begin()) << "\t" << *ic << "\t" << *iw << endl;
      chi2wf.close();
      
      vector <double> weights_GK = mcweights(chi2, npoints, true);
      ofstream chi2wf2((outdirname + "/pdf_GKweights.dat").c_str());

      ic = chi2.begin();
      iw = weights_GK.begin();

      chi2wf2 << "LHAPDF set=   " << lhapdfsetname<<endl;
      chi2wf2 << "Reweight method=   GIELE-KELLER"<<endl;
      chi2wf2 << "ndata=   " << npoints <<endl;
      chi2wf2 << totmc << endl;

      for (; ic != chi2.end(); ic++, iw++)
	chi2wf2 << (ic - chi2.begin()) << "\t" << *ic << "\t" << *iw << endl;
      chi2wf2.close();
      
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
		msg = (string)"S: Error: inconsistent number of MC replica per point";
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
  for (int j = 0; j < totasym; j++)
    {
      for (int i = 0; i < npoints; i++)
	{
	  //account for the sign flip due to applying a theory variation as a shift to the data
	  systasym_.betaasym_[i][0][nsysloc] = -(pointsmap[i].th_asym_p[j] - pointsmap[i].thc) / pointsmap[i].thc;
	  systasym_.betaasym_[i][1][nsysloc] = -(pointsmap[i].th_asym_m[j] - pointsmap[i].thc) / pointsmap[i].thc;
	  systema_.beta_[i][nsysloc] = 0.5*(systasym_.betaasym_[i][0][nsysloc] - systasym_.betaasym_[i][1][nsysloc]);
	  systasym_.omega_[i][nsysloc] = 0.5*(systasym_.betaasym_[i][0][nsysloc] + systasym_.betaasym_[i][1][nsysloc]);
	  
	  //systema_.beta_[i][nsysloc] = (pointsmap[i].th_asym_p[j] - pointsmap[i].th_asym_m[j]) / 2. / pointsmap[i].thc;
	  //systasym_.betaasym_[i][1][nsysloc] = (pointsmap[i].th_asym_m[j] - pointsmap[i].thc) / pointsmap[i].thc;
	  //systasym_.betaasym_[i][0][nsysloc] = (pointsmap[i].th_asym_p[j] - pointsmap[i].thc) / pointsmap[i].thc;
	  //systasym_.omega_[i][nsysloc] = (pointsmap[i].th_asym_p[j] + pointsmap[i].th_asym_m[j] - 2*pointsmap[i].thc) / 2. / pointsmap[i].thc;
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
	  systema_.beta_[i][nsysloc] = - (pointsmap[i].th_hess_s[j] - pointsmap[i].thc) / pointsmap[i].thc;  // negative sign here, opposite to data
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
	  
	  //systema_.beta_[i][nsysloc] = (pointsmap[i].th_env_p - pointsmap[i].th_env_m) / 2.;
	  //systasym_.betaasym_[i][1][nsysloc] = pointsmap[i].th_env_m;
	  //systasym_.betaasym_[i][0][nsysloc] = pointsmap[i].th_env_p;
	  //systasym_.omega_[i][nsysloc] = (pointsmap[i].th_env_p + pointsmap[i].th_env_m) / 2.;
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

  //save the number of total systematic uncertainties including PDF nuisance parameters
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

  // Set PDF nuisance parameter name
  for (int j = systema_.nsys_; j < nsyspdf; j++)
    {
      char nuispar[64];
      sprintf (nuispar, "PDF_nuisance_param_%02d", j+1 - systema_.nsys_);
      int len = strlen(nuispar);
      memset (nuispar+len,' ',64-len);
      strcpy(systema_.system_[j],nuispar);
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
  if (totmc || tothess_s)
    symm = true;
  writetheoryfiles_(nsysloc-systema_.nsys_, theo_cent, symm);

  //Calculate total theory error to be stored in fittedresults.txt for plots
  for (map <int, point>::iterator  pit = pointsmap.begin(); pit != pointsmap.end(); pit++)
    {
      pit->second.th_err_up = 0;
      pit->second.th_err_dn = 0;
    }
  if (totmc)
    for (map <int, point>::iterator pit = pointsmap.begin(); pit != pointsmap.end(); pit++)
      {
	pit->second.th_err_up = pit->second.th_mc_var;
	pit->second.th_err_dn = pit->second.th_mc_var;
      }
  if (totasym)
    for (int i = 0; i < npoints; i++)
      {
	vector <double> xi;
	xi.push_back(pointsmap[i].thc);
	for (int j = 0; j < totasym; j++)
	  {
	    xi.push_back(pointsmap[i].th_asym_p[j]);
	    xi.push_back(pointsmap[i].th_asym_m[j]);
	  }
	double eplus, eminus;
	ahessdeltaasym(xi, eplus, eminus);
	pointsmap[i].th_err_up = eplus;
	pointsmap[i].th_err_dn = eminus;
      }
  if (tothess_s > 0)
    for (int i = 0; i < npoints; i++)
      {
	vector <double> xi;
	xi.push_back(pointsmap[i].thc);
	for (int j = 0; j < tothess_s; j++)
	  xi.push_back(pointsmap[i].th_hess_s[j]);
	pointsmap[i].th_err_up = shessdelta(xi);
	pointsmap[i].th_err_dn = shessdelta(xi);
      }
  //Square add parametrisation envelope
  if (totpar > 0)
    for (int i = 0; i < npoints; i++)
      {
	pointsmap[i].th_err_up = sqrt(pow(pointsmap[i].th_err_up,2) + pow(pointsmap[i].th_env_p-pointsmap[i].thc,2));
	pointsmap[i].th_err_dn = sqrt(pow(pointsmap[i].th_err_dn,2) + pow(pointsmap[i].th_env_m-pointsmap[i].thc,2));
      }

  //Square add scale variations
  if (totscale)
    for (int i = 0; i < npoints; i++)
      {
	vector <double> xi;
	xi.push_back(pointsmap[i].thc);
	for (int j = 0; j < totscale; j++)
	  {
	    xi.push_back(pointsmap[i].th_scale_p[j]);
	    xi.push_back(pointsmap[i].th_scale_m[j]);
	  }
	double eplus, eminus;
	ahessdeltaasym(xi, eplus, eminus);
	pointsmap[i].th_err_up = sqrt(pow(pointsmap[i].th_err_up,2) + pow(eplus,2));
	pointsmap[i].th_err_dn = sqrt(pow(pointsmap[i].th_err_up,2) + pow(eminus,2));
      }

  //Set Error in fortran common block
  for (map <int, point>::iterator pit = pointsmap.begin(); pit != pointsmap.end(); pit++)
    {
      c_theo_.theo_tot_up_[pit->first] = pit->second.th_err_up;
      c_theo_.theo_tot_down_[pit->first] = pit->second.th_err_dn;
      if (clhapdf_.scale68_)
	{
	  c_theo_.theo_tot_up_[pit->first] /= 1.64;
	  c_theo_.theo_tot_down_[pit->first] /= 1.64;
	}
    }
  
  if (pdfprofile || scaleprofile)
    {
      cout << "-------------------------------------------" << endl;
      cout << "Chi2 test on central prediction with PDF and/or scale uncertainties:" << endl;

      // Add PDF nuisance parameters
      for (int i = systema_.nsys_; i < nsyspdf; i++)
	{
	  sysmeas_.n_syst_meas_[i] = npoints; //PDF systematic uncertainties apply to all points
	  for (int j = 0; j < npoints; j++)
	    sysmeas_.syst_meas_idx_[i][j] = j + 1;
	  systscal_.sysscalingtype_[i] = 1;  //Apply linear scaling to PDF uncertainties
	  //      systscal_.sysscalingtype_[i] = 0;  //No scaling for PDF uncertainties
	  if ((nsyspdf-i) <= clhapdf_.nremovepriors_)
	    {
	      char nuispar[64];
	      sprintf (nuispar, "PDF_nuisance_param_%02d", i+1 - systema_.nsys_);
	      cout << "Remove prior for syst " << nuispar << endl;
	      string msg = (string) "I: Remove prior for systematic " + nuispar;
	      hf_errlog_(15082401, msg.c_str(), msg.size());
	      csystprior_.syspriorscale_[i] = 0.;
	    }
	  csysttype_.isysttype_[i] = 2; // THEORY
	}

      // Add QCD scale nuisance parameters
      for (int i = nsyspdf; i < nsysloc; i++)
	{
	  sysmeas_.n_syst_meas_[i] = npoints; //QCD scale systematic uncertainties apply to all points
	  for (int j = 0; j < npoints; j++)
	    sysmeas_.syst_meas_idx_[i][j] = j + 1;
	  systscal_.sysscalingtype_[i] = 1;  //Apply linear scaling to QCD scale uncertainties
	  //      systscal_.sysscalingtype_[i] = 0;  //No scaling
	  //To remove the prior on the scale uncertainties parameters, uncomment the following lines
	  /*
	  char nuispar[64];
	  sprintf (nuispar, "scale_nuisance_param_%02d", i+1 - nsyspdf);
	  cout << "Remove prior for syst " << nuispar << endl;
	  string msg = (string) "I: Remove prior for systematic " + nuispar;
	  hf_errlog_(15082401, msg.c_str(), msg.size());
	  csystprior_.syspriorscale_[i] = 0.;
	  */
	  csysttype_.isysttype_[i] = 2; // THEORY
	}

      //Add the PDF and/or scale uncertainties to the total number of systematic uncertainties
      systema_.nsys_ = nsysloc;
      systematicsflags_.resetcommonsyst_ = true;

      cout << "Total Number of systematic uncertainties: " << systema_.nsys_ << endl;
      fname = outdirname + "/Results.txt";
      fopen_(85, fname.c_str(), fname.size());
      chi2tot = chi2data_theory_(3);
      fclose_(85);
    }
  else
    {
      bool cp = system(((string)"cp " + outdirname + "/fittedresults.txt " 
			+ outdirname + "/fittedresults.txt_set" + tag).c_str());
    }
}
#endif
