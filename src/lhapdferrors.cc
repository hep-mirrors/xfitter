#include "herafitter_cpp.h"

#ifdef LHAPDF_ENABLED
#include <LHAPDF/LHAPDF.h>
#endif

#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <math.h>


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
  double th_mc_sum2;         //monte carlo square sum
  double th_mc_var;          //monte carlo var
  double th_err_up;          //theory error up
  double th_err_dn;          //theory error down
};

//return error if LHAPDF is not enabled
#ifndef LHAPDF_ENABLED
void get_lhapdferrors_()
{
  string msg = "S: Call to LHAPDFErrors but LHAPDF is not enabled. Run ./configure --enable-lhapdf and link the executable";
  hf_errlog_(14060101, msg.c_str(), msg.size());
}
#else

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
  int MonteCarloPDFErr = 0;
  int AsymHessPDFErr = 0;
  int SymmHessPDFErr = 0;
  int ModPDFErr = 0;
  int ParPDFErr = 0;

  map <int, point> pointsmap;

  string lhapdfset = string(clhapdf_.lhapdfset_, 128);
  lhapdfset = lhapdfset.erase(lhapdfset.find_last_not_of(" ")+1, string::npos);
  string lhapdfvarset = string(clhapdf_.lhapdfvarset_, 128);
  lhapdfvarset = lhapdfvarset.erase(lhapdfvarset.find_last_not_of(" ")+1, string::npos);

  string outdirname = string(coutdirname_.outdirname_, 128);
  outdirname = outdirname.erase(outdirname.find_last_not_of(" ")+1, string::npos);

  getpdfunctype_heraf_(lhapdfset.c_str(), MonteCarloPDFErr, AsymHessPDFErr, SymmHessPDFErr);
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

  //Main loop on PDF members, and on EIG and VAR pdf if available
  int cset = 0; //counter on total pdf members;
  int isys = 1; //counter on PDF error index;
  for (int pdfset = 0; pdfset < 2; pdfset++)
    {
      if (pdfset ==  0)
	{
	  LHAPDF::initPDFSet(lhapdfset.c_str());
	  //	  LHAPDF::initPDFSetByName(lhapdfset.c_str());
	  getpdfunctype_heraf_(lhapdfset.c_str(), MonteCarloPDFErr, AsymHessPDFErr, SymmHessPDFErr);
	}
      else if (pdfset ==  1)
	{
	  if (lhapdfvarset == "")
	    continue;
	  LHAPDF::initPDFSet(lhapdfvarset.c_str());
	}
      //Number of PDF members
      int nsets = LHAPDF::numberPDF(); 	  
      cout << "Number of PDF members for this set: " << nsets << endl;

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
	      if (iset > (nsets - clhapdf_.nparvar_))
		ParPDFErr = true;
	      else
		ModPDFErr = true;
	    }

	  //Evaluate chi2 for the current PDF member
	  double chi2tot;
	  int iflag = min(2,max(1, cset+1));
	  chi2tot = chi2data_theory_(iflag);
	  char chi2c[20];
	  sprintf(chi2c, "%.2f", chi2tot);
	  cout << setw(20) << "PDF set number: " << setw(5) << iset 
	       << setw(15) << "chi2=" << chi2c 
	       << setw(15) << "ndf=" << cfcn_.ndfmini_ 
	       << endl;

	  //Write results of chi2 calculation
	  char tag[10];
	  sprintf (tag, "_%04d", cset);
	  writefittedpoints_(); //write out fittedresults.txt file
	  bool mv = system(((string)"mv " + outdirname + "/fittedresults.txt " 
			    + outdirname + "/fittedresults.txt_set" + tag).c_str());

	  //Store PDF members
	  string filename = outdirname + "/pdfs_q2val_";
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
		  filename += "m_";  //Asymmetric minus suffix
		else
		  filename += "p_";  //Asymmetric plus suffix
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

	  //Store all PDF variations
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
		    pointsmap[i].th_asym_m.push_back(c_theo_.theo_[i]);
		  else
		    pointsmap[i].th_asym_p.push_back(c_theo_.theo_[i]);
		else if (ModPDFErr)
		  if ((cset%2) == 1)
		    pointsmap[i].th_asym_m.push_back(c_theo_.theo_[i]);
		  else
		    pointsmap[i].th_asym_p.push_back(c_theo_.theo_[i]);
		else if (ParPDFErr)
		  pointsmap[i].th_par.push_back(c_theo_.theo_[i]);
	      }
	  
	  if (iset != 0)
	    if (MonteCarloPDFErr || SymmHessPDFErr || ParPDFErr)
	      isys++;
	    else if ((iset%2) == 0) //set the same index for Up and Down variation of asymmetric PDF errors
	      isys++;
	  cset++;
	} //End of loop on PDF members
    } //End loop on eig var PDF sets


  //Monte Carlo replica PDF uncertainties
  int totmc = 0;
  for (map <int, point>::iterator  pit = pointsmap.begin(); pit != pointsmap.end(); pit++)
    {
      //MC replica mean
      pit->second.th_mc_mean = 0;
      pit->second.th_mc_sum2 = 0;
      for (vector <double>::iterator mcit = pit->second.th_mc.begin(); mcit != pit->second.th_mc.end(); mcit++)
	{
	  pit->second.th_mc_mean += *mcit;
	  pit->second.th_mc_sum2 += pow(*mcit, 2);
	}
      pit->second.th_mc_mean /= pit->second.th_mc.size();
      pit->second.th_mc_var = sqrt(pit->second.th_mc_sum2/pit->second.th_mc.size() 
				   - pow(pit->second.th_mc_mean, 2)) / pit->second.th_mc_mean;
      totmc = pit->second.th_mc.size();
    }
  if (totmc > 0)
    {
      char num[10];
      sprintf (num, "%d", totmc);
      msg = (string) "I: Found " + num + " Monte Carlo PDF uncertainties variations";
      hf_errlog_(25051401, msg.c_str(), msg.size());
    }
  //Evaluate MC covariance matrix
  int dim = npoints;
  double covmx[dim*dim];
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
  if (totmc > 0)
    {
      double beta_from_covmx[dim*npoints];
      double alpha_from_covmx[dim];
      int ncorr = 0;	
      getnuisancefromcovar_(dim,npoints,npoints,
			    covmx,beta_from_covmx,0,
			    ncorr,alpha_from_covmx);

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
    }

  //Asymmetric PDF uncertainties, including hessian and model
  int totasym = 0;
  for (map <int, point>::iterator  pit = pointsmap.begin(); pit != pointsmap.end(); pit++)
    if (pit->second.th_asym_p.size() != pit->second.th_asym_m.size())
      {
	msg = (string)"S: Error: inconsistent number of positive and negative asymmetric PDF variations, check your NPARVAR setting";
	hf_errlog_(25051403, msg.c_str(), msg.size());
      }
    else
      totasym = pit->second.th_asym_p.size(); //include both model and parametrisation uncertainties
  if (totasym > 0)
    {
      char num[10];
      sprintf (num, "%d", totasym);
      msg = (string) "I: Found " + num + " asymmetric PDF uncertainties variations";
      hf_errlog_(25051404, msg.c_str(), msg.size());
    }
  for (int j = 0; j < totasym; j++)
    {
      for (int i = 0; i < npoints; i++)
	{
	  systema_.beta_[i][nsysloc] = (pointsmap[i].th_asym_p[j] - pointsmap[i].th_asym_m[j]) / 2. / pointsmap[i].thc;
	  systasym_.betaasym_[i][1][nsysloc] = (pointsmap[i].th_asym_m[j] - pointsmap[i].thc) / pointsmap[i].thc;
	  systasym_.betaasym_[i][0][nsysloc] = (pointsmap[i].th_asym_p[j] - pointsmap[i].thc) / pointsmap[i].thc;
	  systasym_.omega_[i][nsysloc] = (pointsmap[i].th_asym_p[j] + pointsmap[i].th_asym_m[j] - 2*pointsmap[i].thc) / 2. / pointsmap[i].thc;
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
      msg = (string) "I: Found " + num + " Symmetric hessian PDF uncertainties variations";
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
      msg = (string) "I: Found " + num + " Parametrisation uncertainties";
      hf_errlog_(25051406, msg.c_str(), msg.size());
    }
  //make envelope;
  for (map <int, point>::iterator  pit = pointsmap.begin(); pit != pointsmap.end(); pit++)
    {
      pit->second.th_env_p = 0;
      pit->second.th_env_m = 0;
      for (vector <double>::iterator it = pit->second.th_par.begin(); it != pit->second.th_par.end(); it++)
	{
	  pit->second.th_env_p = max(0.,max((*it - pit->second.thc)/ pit->second.thc, pit->second.th_env_p));
	  pit->second.th_env_m = min(0.,min((*it - pit->second.thc)/ pit->second.thc, pit->second.th_env_m));
	}
    }
  if (totpar > 0)
    {
      for (int i = 0; i < npoints; i++)
	{
	  systema_.beta_[i][nsysloc] = (pointsmap[i].th_env_p - pointsmap[i].th_env_m) / 2.;
	  systasym_.betaasym_[i][1][nsysloc] = pointsmap[i].th_env_m;
	  systasym_.betaasym_[i][0][nsysloc] = pointsmap[i].th_env_p;
	  systasym_.omega_[i][nsysloc] = (pointsmap[i].th_env_p + pointsmap[i].th_env_m) / 2.;
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

  // Set nuisance parameter name
  for (int j = systema_.nsys_; j < nsysloc; j++)
    {
      char nuispar[64];
      sprintf (nuispar, "PDF_nuisance_param_%02d", j+1 - systema_.nsys_);
      int len = strlen(nuispar);
      memset (nuispar+len,' ',64-len);
      strcpy(systema_.system_[j],nuispar);
    }

  //Set central theory value in fortran common block
  double theo_cent[2500];
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
	double err2_up = 0;
	double err2_dn = 0;
	for (int j = 0; j < totasym; j++)
	  {
	    err2_up += pow(max(max((pointsmap[i].th_asym_p[j] - pointsmap[i].thc) / pointsmap[i].thc,
				   (pointsmap[i].th_asym_m[j] - pointsmap[i].thc) / pointsmap[i].thc), 0.),2);
	    err2_dn += pow(max(max((pointsmap[i].thc - pointsmap[i].th_asym_p[j]) / pointsmap[i].thc,
				   (pointsmap[i].thc - pointsmap[i].th_asym_m[j]) / pointsmap[i].thc), 0.),2);
	  }
	pointsmap[i].th_err_up = sqrt(err2_up);
	pointsmap[i].th_err_dn = sqrt(err2_dn);
      }
  if (tothess_s > 0)
    for (int i = 0; i < npoints; i++)
      {
	double err2_up = 0;
	double err2_dn = 0;
	for (int j = 0; j < tothess_s; j++)
	  {
	    err2_up += pow((pointsmap[i].th_hess_s[j] - pointsmap[i].thc) / pointsmap[i].thc, 2);
	    err2_dn += pow((pointsmap[i].th_hess_s[j] - pointsmap[i].thc) / pointsmap[i].thc, 2);
	  }
	pointsmap[i].th_err_up = sqrt(err2_up);
	pointsmap[i].th_err_dn = sqrt(err2_dn);
      }
  //Square add parametrisation envelope
  if (totpar > 0)
    for (int i = 0; i < npoints; i++)
      {
	pointsmap[i].th_err_up = sqrt(pow(pointsmap[i].th_err_up,2) + pow(pointsmap[i].th_env_p,2));
	pointsmap[i].th_err_dn = sqrt(pow(pointsmap[i].th_err_dn,2) + pow(pointsmap[i].th_env_m,2));
      }

  //Set Error in fortran common block
  for (map <int, point>::iterator pit = pointsmap.begin(); pit != pointsmap.end(); pit++)
    {
      c_theo_.theo_tot_up_[pit->first] = pit->second.th_err_up * pit->second.thc;
      c_theo_.theo_tot_down_[pit->first] = pit->second.th_err_dn * pit->second.thc;
      if (clhapdf_.scale68_)
	{
	  c_theo_.theo_tot_up_[pit->first] /= 1.64;
	  c_theo_.theo_tot_down_[pit->first] /= 1.64;
	}
    }
  
  LHAPDF::initPDF(0);
  c_alphas_.alphas_ = LHAPDF::alphasPDF(boson_masses_.mz_);

  cout << "-------------------------------------------" << endl;
  cout << "Chi2 test on central prediction:" << endl;
  string fname = outdirname + "/Results_00.txt";
  fopen_(85, fname.c_str(), fname.size());
  double chi2tot = chi2data_theory_(3);
  fclose_(85);

  cout << "-------------------------------------------" << endl;
  cout << "Chi2 test on central prediction with PDF uncertainties:" << endl;

  for (int i = systema_.nsys_; i < nsysloc; i++)
    {
      sysmeas_.n_syst_meas_[i] = npoints; //PDF systematic uncertainties apply to all points
      for (int j = 0; j < npoints; j++)
	sysmeas_.syst_meas_idx_[i][j] = j + 1;
      systscal_.sysscalingtype_[i] = 1;  //Apply linear scaling to PDF uncertainties
    }

  //Add the PDF uncertainties to the total number of systematic uncertainties
  systema_.nsys_ = nsysloc;
  systematicsflags_.resetcommonsyst_ = true;

  cout << "Total Number of systematic uncertainties: " << systema_.nsys_;
  fname = outdirname + "/Results.txt";
  fopen_(85, fname.c_str(), fname.size());
  chi2tot = chi2data_theory_(3);
  fclose_(85);
}
#endif
