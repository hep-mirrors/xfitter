#include <map>
#include <iostream>
#include <sstream>
#include <string>
using std::cout;
using std::endl;
using namespace std;

#include <cstdlib>
#include <sys/stat.h>
#include <cmath>

#include "pdfs.h"
#include "utils.h"

// ROOT
#include "TStyle.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1F.h"
#include "TColor.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TRandom3.h"

// LHAPDF
#include "LHAPDF/LHAPDF.h"

using namespace std;


/// externally defined alpha_s and pdf routines for fortran 
/// callable convolution wrapper
extern "C" void nnnpdf_reweight_( int &npdfs, const char* s) {

  int size = strlen(s);

  char *end = (char*)s + size - 1;
  while (end >= s && isspace(*end))
    end--;
  *(end + 1) = '\0';

  while (*s && isspace(*s))
    s++;
 
  size_t NPOINTS=100;

  cout << "READING IN NNPDF steering FILE: " << s << endl; 

  // ANALYSIS PARAMETERS
  // Reweighting parameters
  rwparam rpar;
  parse_param_input(s,rpar); //------------> chi2 is read in...
  
  // Prior PDF parameters
  const string NAME = rpar.prior;
    
  // PDF parameters
  PDFparams par;
  par.nrep = npdfs;
  par.nflav = 13;
  par.nx=0;
  par.NAME=NAME;
	
  PDF xf(par);
	
  // ************************FILE HANDLING **********************************

  ostringstream filename;
           
  // Weight Histogram plotting 
  filename.str("");
  filename <<rpar.outdir<< "/whist-"<<rpar.outdesc<<"."<<rpar.plotform;
  string whistfile=filename.str();

  filename.str("");
  filename <<rpar.outdir<< "/palpha-"<<rpar.outdesc<<"."<<rpar.plotform;
  string paout=filename.str();
	
  //****************************WEIGHTS**************************************  
	
  // Reweight		
  xf.Reweight(rpar);
  xf.CheckWeights();	
	
  // Unweight (for LHGrid out)
  PDF uxf(xf,rpar.size);
  // Write out the Reweighted LHGrid
  if (rpar.lhgrid)
    uxf.Export(rpar);	
	
  // *********************PLOTTING***************************************
	
  // Weight Histogram
  TCanvas *wHc = new TCanvas ("wHistC","Weight Histogram",12,38,699,499);
  TH1F* wHist= new TH1F("wH", "Weight Histogram;Weight;Frequency;", 50, -7, 1);
  
  // Log binning  
  BinLogX(wHist);
  wHc->GetPad(0)->SetLogx();
		
  vector<double> w=xf.GetWeights();
		
  for (size_t i=0; i<w.size(); i++)
    if (w[i]!=0)
      wHist->Fill(w[i]);
		
  wHist->Draw();
  wHc->Print(("output/"+whistfile).c_str());
		
  delete wHc;
	
  // P(alpha) plot
  TGaxis::SetMaxDigits(3);
    
  Double_t alpha[NPOINTS],palph[NPOINTS];
  double ptot=0;
  for (size_t i=0; i<NPOINTS; i++)
    {
      alpha[i]=5*(Double_t ) i/(Double_t  )NPOINTS + 0.1; 
      palph[i]=xf.Palpha(alpha[i]);
      ptot=ptot+palph[i];
    }
  
    // Roughly normalise
  double intp= integrate(palph, NPOINTS, 5.0/((double) NPOINTS));
    for (size_t i=0; i<NPOINTS; i++)
        palph[i]=palph[i]/(intp);
  
  TGraph* dpalpha= new TGraph(NPOINTS,alpha,palph);
  TCanvas *dCpa = new TCanvas("paPlot", "P(#alpha)",12,38,699,499);		
    
  //  dCpa->GetPad(0)->SetLogx();
  
  gStyle->SetOptStat(0);
  dCpa->SetBorderSize(0);
  dCpa->SetBorderMode(0);
  dCpa->SetFrameFillColor(0);
  dCpa->SetFrameBorderMode(0);
  dCpa->SetFillColor(0);    
        
  dpalpha->SetMinimum(0.);
  dpalpha->SetTitle("P(#alpha)");
    
  dpalpha -> SetLineColor(kBlue);
  dpalpha -> SetLineWidth(3); 
  dpalpha -> SetLineStyle(7); 
	
  TAxis* xaxpa = dpalpha->GetXaxis();
  xaxpa->SetTitle("#alpha");
  TAxis* yaxpa = dpalpha->GetYaxis();
  yaxpa->SetTitle("P(#alpha)");
    
  yaxpa->SetTitleOffset(1.1);
    
  dpalpha->Draw("AL");
	
  cout << (("output/"+paout).c_str()) << endl;
  cout <<dCpa << endl;
  dCpa->Print(("output/"+paout).c_str());
  cout << "hello"<<endl;
  //  exit(0);
  return;
}


/// callable convolution wrapper: create PDF replicas
extern "C" void create_randompdfreplicas_( const char* pdfset, const char* directory, int * nrep, int  * isSymmetric) {

  int nRandomReplicas = *nrep;
  const size_t cnRandomReplicas = nRandomReplicas+1;
  stringstream ofilename;
  ofilename.str("");
  ofilename << "output/" << directory << "/" << pdfset << Form("_%dInputReplicas",nRandomReplicas)<<"/"<< pdfset << Form("_%dInputReplicas",nRandomReplicas);
  std::cout << "nnpdfrw isSymmetric  " << *isSymmetric<< std::endl;
  int hasSymmErrorSets = *isSymmetric;

  ifstream file((ofilename.str()+"_0000.dat").c_str());
  if(file){
    std::cout<<"\033[1;31mOutput file " <<((ofilename.str()+"_0000.dat").c_str())<<" already exists. Skipping generation of random replicas!\033[0m"<<std::endl;
    return;
  }
  int size = strlen(directory);

  char *end = (char*)directory + size - 1;
  while (end >= directory && isspace(*end))
    end--;
  *(end + 1) = '\0';
  
  while (*directory && isspace(*directory))
    directory++;
  
  size = strlen(pdfset);

  char *end2 = (char*)pdfset + size - 1;
  while (end2 >= pdfset && isspace(*end2))
    end2--;
  *(end2 + 1) = '\0';
  
  while (*pdfset && isspace(*pdfset))
    pdfset++;


  cout << "**********"<< directory<< endl;
  cout << "**********"<< pdfset<< endl;

  const size_t npx=100, npq2=50;
  typedef double xdim[npx][npq2][13];
  
 const LHAPDF::PDFSet set(pdfset);
 /// create also VAR version!
 /*
 bool var=false;
 string setname(pdfset); 
 if (setname.find("EIG")!=std::string::npos) {
  setname.replace(setname.length()-3,3,"VAR");
  cout << "Automatically also testing for set: ** "<< setname<< "  -- in case of crash: Check whether the PDF is available" << endl;
  LHAPDF::PDFSet * varset = new LHAPDF::PDFSet(setname.c_str());
  }*/

  // Init x and q2 grid points
  vector<double> xg, qg;
  double qmin=sqrt(LHAPDF::getQ2min(0));
  double qmax=sqrt(LHAPDF::getQ2max(0));
  double xmin=LHAPDF::getXmin(0);
  double xmax=LHAPDF::getXmax(0);
  
  // FORTRAN STANDARDS
  double XCH=0.1;
  
    // Set up x, q2 grids
  for (size_t i=1; i<(npx+1); i++)
    if (i<50)
      xg.push_back(xmin*pow(XCH/xmin,2.0*(((double) i ) -1)/(((double) npx ) -1)));
    else
      xg.push_back(XCH+(xmax-XCH)*(((double) i ) -51)/(((double) npx ) -51) );
  
  for (size_t i=1; i<(npq2+1); i++)
    qg.push_back(sqrt((qmin*qmin*pow((qmax*qmax)/(qmin*qmin),(((double) i ) -1)/(((double) npq2 ) -1)))));
  
  const int nflavs=13;
  int flavs[nflavs]={-6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 21};
  
  double oas = LHAPDF::getOrderAlphaS();
  double opdf = LHAPDF::getOrderPDF();
  
  cout <<endl<<"Writing out replica LHAPDF grid: "<< pdfset << "_InputReplicas (.LHgrid)" <<endl;
  if (hasSymmErrorSets) cout <<endl<<"creating "<< nRandomReplicas << " random replicas PDF with symmetric errors" <<endl;
  else  cout <<endl<<"creating "<< nRandomReplicas << " random replicas from PDF with asymmetric errors" <<endl;
  cout <<"Using LHAPDF version: "<< LHAPDF::getVersion()<<endl;

  cout << "Creating random replica input LHgrid file" << endl;
  cout << ofilename.str() << endl;
  cout << "------------------------------------------" << endl;

  string dir = string(directory);
   mkdir( ("output/"+dir).c_str() ,0777);
  dir=dir+"/";
  dir=dir+pdfset;
  dir=dir+Form("_%d",nRandomReplicas);
  dir=dir+"InputReplicas";
  mkdir( ("output/"+dir).c_str() ,0777);

  // Create out files: 
  // header: contains header, info on PDF sets etc. --> .info
  // centralval: contains set 0 / mean --> _0000.dat
  // all other sets: _0001.dat - 000X.dat --> later... (single files)

   // Create XXX.info files: contains header, info on PDF sets etc.

  // string outdir = "output/";
  // outdir=outdir+dir;
  ofstream lhaout_header( (ofilename.str()+".info").c_str() );
  //  ofstream lhaout_centralval( (ofilename+"_0000.dat").c_str() );
  // ofstream lhaout_errorsets( (ofilename+"/errorsets.txt").c_str() );

  // Write out LHAPDF preamble
  lhaout_header.precision(7);
  lhaout_header << scientific;

  lhaout_header <<"SetDesc: \"";
  lhaout_header << pdfset<< "  (lhapdf set name) created by the HeraFitter package using " << nRandomReplicas<< " random replicas, mem=0 => average on replicas; mem=1-" << nRandomReplicas << "=> PDF replicas\""<<endl;
  lhaout_header <<"Authors: ... " << endl;
  lhaout_header <<"Reference: ..." << endl;
  lhaout_header <<"Format: lhagrid1" << endl;
  lhaout_header <<"DataVersion: 6" << endl;
  lhaout_header <<"NumMembers: " << (int)(nRandomReplicas+1) << endl;
  lhaout_header <<"Particle: 2212" << endl;
  lhaout_header <<"Flavors: [-6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 21]" << endl;

  if (oas==0.0)
    lhaout_header<< "OrderQCD: 0"<<endl;
  else if (oas==1.0)
    lhaout_header<< "OrderQCD: 1"<<endl;
  else if (oas==2.0)
    lhaout_header<< "OrderQCD: 2"<<endl;
  else
    {
      cout <<"ERR: invalid asorder"<<endl;
      cout <<oas<<endl;
      exit(1);
    }

  lhaout_header <<"FlavorScheme: variable" << endl;
  lhaout_header <<"NumFlavors: 6" << endl;
  lhaout_header <<"ErrorType: replicas" << endl;
  lhaout_header <<"XMin: " << LHAPDF::getXmin(0)<<endl;
  lhaout_header <<"XMax: " <<LHAPDF::getXmax(0)<< endl;
  lhaout_header <<"QMin: " <<sqrt(LHAPDF::getQ2min(0))<< endl;
  lhaout_header <<"QMax: " << sqrt(LHAPDF::getQ2max(0))<<endl;
 lhaout_header <<"MZ: 91.1876" << endl;
 lhaout_header <<"MUp: 0" << endl;
 lhaout_header <<"MDown: 0" << endl;
 lhaout_header <<"MStrange: 0" << endl;
 lhaout_header <<"MCharm: " << LHAPDF::getQMass(4)<<endl;
 lhaout_header <<"MBottom: "<< LHAPDF::getQMass(5)<<endl;
 lhaout_header <<"MTop: "<< LHAPDF::getQMass(6)<<endl;
 lhaout_header <<"AlphaS_MZ: " <<LHAPDF::alphasPDF(91.1876)<< endl;

  if (oas==0.0)
    lhaout_header<< "AlphaS_OrderQCD: 0"<<endl;
  else if (oas==1.0)
    lhaout_header<< "AlphaS_OrderQCD: 1"<<endl;
  else if (oas==2.0)
    lhaout_header<< "AlphaS_OrderQCD: 2"<<endl;
  else
    {
      cout <<"ERR: invalid asorder"<<endl;
      cout <<oas<<endl;
      exit(1);
    }

   lhaout_header <<"AlphaS_Type: ipol" << endl;
  lhaout_header <<"AlphaS_Qs: [" ;
  for (size_t i=0; i<npq2; i++) {
    lhaout_header<<" "<<  qg[i] ;// qg --> qgrid --> points in q
  }
  lhaout_header <<"]" << endl;     

 lhaout_header <<"AlphaS_Vals: [" ;
  for (size_t i=0; i<npq2; i++) {
    lhaout_header<<" "<<  LHAPDF::alphasPDF(qg[i]) ; // qg --> qgrid --> points in q
  }
  lhaout_header <<"]" << endl;     
 
  lhaout_header << "AlphaS_Lambda4:  "<<LHAPDF::getLam4(0)<<endl;
  lhaout_header << "AlphaS_Lambda5:  "<<LHAPDF::getLam5(0)<<endl;
  lhaout_header.close();

  // Now write out single error set files
  // first: preparation of arrays / random numbers
  // Determine values of xfx, place in xfxval array, compute average over replicas for 0th member PDF
  vector<double> xfxavg;
  const size_t nxdim=2;
  xdim *xfxval = new xdim[nxdim];

  TRandom3 * random=new TRandom3() ;

  double *rand = new double[(LHAPDF::numberPDF()+1)*nRandomReplicas+1]; // array of random numbers
  for (int ir=0;ir<=(LHAPDF::numberPDF()+1)*nRandomReplicas;ir++) { // 20*40+1 random numbers 
    rand[ir] = random->Gaus(0.0,1.0); // Gaussian random number
  }
 
  const vector<LHAPDF::PDF*> pdfs = set.mkPDFs();
 
  cout << "USING NREP: "<<nRandomReplicas<<endl;
  // Write out the contents of the xfxvals array to the LHgrid
  for (int n=0; n<nRandomReplicas; n++)
  {
    cout <<"\n\033[F\033[J";
   cout << "Writing replica: "<<n+1<<"/"<<nRandomReplicas<<endl;
    cout<<"\n\033[F\033[J";

    string filenumber="";
    if (n+1<10) filenumber="000";
    else if (n+1<100) filenumber="00";
    else if (n+1<1000) filenumber="0";

    cout <<  ofilename.str()<<"_"<<filenumber<<(n+1)<<".dat"<< endl;
    string errset_filename=ofilename.str()+"_";
    errset_filename=errset_filename+filenumber;
    errset_filename=errset_filename+Form("%d",n+1);
    ofstream lhaout_errorset( (errset_filename+".dat").c_str() );
    lhaout_errorset<< "PdfType: replica"<<endl;
    lhaout_errorset<< "Format: lhagrid1"<<endl;
    lhaout_errorset<< "---"<<endl;

 // rest in fixed 
  lhaout_errorset.precision(8);
  lhaout_errorset << fixed;
  lhaout_errorset << scientific;
 
    /// x grid
    /// q grid
    /// flavour array

 for (size_t i=0; i<npx; i++) {
   lhaout_errorset<<" "<<  xg[i]; // xg --> xgrid --> points in x   lhaout_header<<  " "<<xg[i]<<endl; 
  }
 lhaout_errorset<< endl;

 for (size_t i=0; i<npq2; i++) {
   lhaout_errorset<<" "<<  qg[i]; // qg --> qgrid --> points in q
  }
 lhaout_errorset<< endl;

lhaout_errorset<< "-6 -5 -4 -3 -2 -1 1 2 3 4 5 6 21" << endl;
 
    for (size_t x=0; x<npx; x++)
    {
      //      if ((x+1)%10==0) {
      cout << Form("\r \033[1;34mGenerating x-Q2 grid for replica %d/%d -- %d",n+1,nRandomReplicas, (int)(100*(x+1.0)/npx) )<<" % done \033[0m" << std::flush; //<< endl;

      for (size_t q=0; q<npq2; q++)
	    {
	      double xval=xg[x], qval=qg[q];
	      for (int i=0; i<nflavs; i++) 
        {
	  //          pdfs[imem]->xfxQ(21,x,Q); 
          double centralvalue_FS0 = pdfs[0]->xfxQ(flavs[i],xval,(qval));//LHAPDF::xfx(xval,(qval),flavs[i]);//
          double val = centralvalue_FS0;
	     if ( ! (hasSymmErrorSets) ) {  // loop over actual e eigenvectors in order to create random replica n from them

            for (int e=1; e<LHAPDF::numberPDF()+1; e++)       {              // asymmetric eigenvector sets
              if (e%2==1) {		  
                double r = rand[(n)*LHAPDF::numberPDF()+e];
		double eigenvector_pos=pdfs[e]->xfxQ(flavs[i],xval,(qval));//LHAPDF::xfx(xval,(qval),flavs[i]);//LHAPDF::xfx(i,xval,(qval));
		double eigenvector_neg=pdfs[e+1]->xfxQ(flavs[i],xval,(qval));//LLHAPDF::xfx(xval,(qval),flavs[i]);//LHAPDF::xfx(i,xval,(qval)) ;
		val=val + 0.5*(eigenvector_pos - eigenvector_neg)*(r);
              }        }  	     }
          else {
            for (int e=1; e<LHAPDF::numberPDF()+1; e++) { // symmetric eigenvector sets
	      double r = rand[(n)*LHAPDF::numberPDF()+e];
              LHAPDF::initPDF( e );
              double eigenvector=LHAPDF::xfx(xval,(qval),flavs[i]);//LHAPDF::xfx(i,xval,(qval));
              val=val+ (eigenvector - centralvalue_FS0) * (r); 
           }	  }

          if (abs(val)<1e-10) val=0;

 	  xfxval[0][x][q][(i)]=val;
	   
          xfxval[1][x][q][(i)]=xfxval[1][x][q][(i)]+val;

          lhaout_errorset << " "<<xfxval[0][x][q][(i)]<<"  "; 
        }
	      lhaout_errorset <<endl; 
           
      }
    }
    lhaout_errorset<< "---"<<endl;
    lhaout_errorset.close();

  }///for (int n=0; n<nRandomReplicas; n++) 

  // now create average file

  ofstream lhaout_centralval( (ofilename.str()+"_0000.dat").c_str() );

   lhaout_centralval<< "PdfType: central"<<endl;
    lhaout_centralval<< "Format: lhagrid1"<<endl;
    lhaout_centralval<< "---"<<endl;

 // rest in fixed 
  lhaout_centralval.precision(7);
  lhaout_centralval<< fixed;
  lhaout_centralval<< scientific;

    /// x grid
    /// q grid
    /// flavour array
 
 for (size_t i=0; i<npx; i++) {
   lhaout_centralval<<" "<<  xg[i]; // xg --> xgrid --> points in x   lhaout_header<<  " "<<xg[i]<<endl; 
  }
 lhaout_centralval<< endl;

for (size_t i=0; i<npq2; i++) {
   lhaout_centralval<<" "<<  qg[i]; // qg --> qgrid --> points in q
  }
 lhaout_centralval<< endl;

lhaout_centralval<< "-6 -5 -4 -3 -2 -1 1 2 3 4 5 6 21" << endl;

  for (size_t x=0; x<npx; x++)
	{
	  for (size_t q=0; q<npq2; q++)
    {
      double xval=xg[x], qval=qg[q];
      for (int i=0; i<nflavs; i++) 
      {
        lhaout_centralval << " "<<xfxval[1][x][q][(i)]/nRandomReplicas<<"  ";  
      }              
		  lhaout_centralval <<endl; 
    }
	}

    lhaout_centralval<< "---"<<endl;
  lhaout_centralval.close();

  cout <<"LHAPDF Writeout successful" <<endl<<endl;
  delete [] xfxval;

  return;
}


  ////////////////////////////////////////////////////////////////
  /* OLD LHAPDF-5 code 

  string outdir = "output/"+dir;
  ofstream lhaout_header( (outdir+"/header.txt").c_str() );
  ofstream lhaout_centralval( (outdir+"/centralval.txt").c_str() );
  ofstream lhaout_errorsets( (outdir+"/errorsets.txt").c_str() );

  // remove("output/RandomReplicaPDF.LHgrid");

  // system( "cd output");
  // cout << Form("ln -s %s/%s_%dInputReplicas.LHgrid RandomReplicaPDF.LHgrid ",directory,pdfset, nRandomReplicas) << endl;
  // system( Form("ln -s %s/%s_%dInputReplicas.LHgrid RandomReplicaPDF.LHgrid ",directory,pdfset, nRandomReplicas)  );
  // system( Form("ls -rtlh output/%s/%s_%dInputReplicas.LHgrid",directory,pdfset, nRandomReplicas)  );
  // system( "cd ../");
 
  // Write out LHAPDF preamble
  lhaout_header.precision(8);
  lhaout_header << scientific;
  
  lhaout_header<<" \'Version\' \'"<<LHAPDF::getVersion()<<"\'"<<endl;
  lhaout_header <<" \'Description:\'"<<endl;
  lhaout_header << " \' " << pdfset<< " " << nRandomReplicas<< " random replicas\'"<<endl;
  lhaout_header << " \' created by the HeraFitter package:\'"<<endl;
  lhaout_header << " \' arvix number\'"<<endl;
  
  lhaout_header<< " \'Alphas:\'"<<endl;
  if (oas==0.0)
    lhaout_header<< " \'Variable', \'lo\', \'EvolCode\'"<<endl;
  else if (oas==1.0)
    lhaout_header<< " \'Variable', \'nlo\', \'EvolCode\'"<<endl;
  else if (oas==2.0)
    lhaout_header<< " \'Variable', \'nnlo\', \'EvolCode\'"<<endl;
  else
    {
      cout <<"ERR: invalid asorder"<<endl;
      cout <<oas<<endl;
      exit(1);
    }
  
  lhaout_header << " 1 ,   91.2      ,   "<<LHAPDF::getQMass(4)<<"      ,   "<<LHAPDF::getQMass(5)<<"      ,   "<<LHAPDF::getQMass(6)<<endl;     
  
  lhaout_header<< " \'MinMax:\'"<<endl;
  lhaout_header <<" "<<  nRandomReplicas  <<",  1"<<endl;
  lhaout_header <<" "<<LHAPDF::getXmin(0)<<" , "<<LHAPDF::getXmax(0)<<" , "<<LHAPDF::getQ2min(0)<<" , "<<LHAPDF::getQ2max(0)<<endl; 
  
  lhaout_header<< " \'QCDparams:\'"<<endl;
  lhaout_header <<" "<<nRandomReplicas<<",  1"<<endl;
  lhaout_header << " "<<LHAPDF::getLam4(0)<<" ,  "<<LHAPDF::getLam5(0)<<endl;
  
  // alphas values for each member
  const double mz = 91.2;
  lhaout_header<< " \'Parameterlist:\'"<<endl;
  lhaout_header<< " \'list', "<<nRandomReplicas<<" , "<<" 1"<<endl;
  for (size_t i=0; i<(cnRandomReplicas); i++)
    lhaout_header << " "<<LHAPDF::alphasPDF(mz)<<endl;
  
  // Order of evolution
  lhaout_header << " \'Evolution:\'"<<endl;
  
  if (opdf==0.0)
    lhaout_header<< " \'lo\', "<<LHAPDF::getQ2min(0)<<" , "<<1<<endl;
  else if (opdf==1.0)
    lhaout_header<< " \'nlo\', "<<LHAPDF::getQ2min(0)<<" , "<<1<<endl;
  else if (opdf==2.0)
    lhaout_header<< " \'nnlo\', "<<LHAPDF::getQ2min(0)<<" , "<<1<<endl;
  else
    {
      cerr <<"ERR: invalid pdf order"<<endl;
      cerr <<opdf<<endl;
      exit(1);
    }
  
  lhaout_header <<  " \'NNPDF20int\'"<<endl; 
  lhaout_header <<  " "<<nRandomReplicas<<", "<<1<<endl;
  
  
  // Write out x, q2 grid values
  lhaout_header <<fixed<< " "<<npx <<endl;
  
  // x in scientific
  lhaout_header.precision(18);
  lhaout_header << scientific;
  
  for (size_t i=0; i<npx; i++)
    lhaout_header<<  " "<<xg[i]<<endl;
  
  lhaout_header <<" "<<npq2 <<endl;
  
  // Q in fixed
  lhaout_header.precision(18);
  lhaout_header << scientific;
  lhaout_header <<" "<<qg[0]<<endl;
  for (size_t i=0; i<npq2; i++)
    lhaout_header<<" "<<  qg[i]<<endl;
  
  lhaout_header <<" "<< nRandomReplicas <<endl;
  
  // rest in fixed 
  lhaout_errorsets.precision(8);
  lhaout_errorsets << fixed;
  lhaout_centralval.precision(8);
  lhaout_centralval << fixed;
  
  // Determine values of xfx, place in xfxval array, compute average over replicas for 0th member PDF
  vector<double> xfxavg;
  const size_t nxdim=2;
  xdim *xfxval = new xdim[nxdim];

  TRandom3 random;

  double *rand = new double[LHAPDF::numberPDF()*nRandomReplicas+1]; // array of random numbers
  for (int ir=0;ir<=LHAPDF::numberPDF()*nRandomReplicas;ir++) { // 20*40+1 random numbers 
    rand[ir] = random.Gaus(0.0,1.0); // Gaussian random number
  }
  
  cout << "USING NREP: "<<nRandomReplicas<<endl;
  // Write out the contents of the xfxvals array to the LHgrid
  for (int n=0; n<nRandomReplicas; n++)
  {
    cout << "Writing replica: "<<n+1<<"/"<<nRandomReplicas<<endl;
    cout<<"\n\033[F\033[J";

    for (size_t x=0; x<npx; x++)
    {
      cout << "Generating grid: "<<x+1<<"/"<<npx;
      cout<<"\n\033[F\033[J";
      for (size_t q=0; q<npq2; q++)
	    {
	      double xval=xg[x], qval=qg[q];
	      for (int i=-6; i<7; i++) 
        {
          LHAPDF::initPDF( 0 ); 
          double centralvalue_FS0 = LHAPDF::xfx(xval,sqrt(qval),i);
          double val = centralvalue_FS0;
	     if ( ! (hasSymmErrorSets) ) {  // loop over actual e eigenvectors in order to create random replica n from them
            for (int e=1; e<LHAPDF::numberPDF()+1; e++)       {              // asymmetric eigenvector sets
              if (e%2==1) {		  
                double r = rand[(n-1)*LHAPDF::numberPDF()+e];
		double  eigenvector=0;
		if (r<0) {
		  LHAPDF::initPDF( e );
		  eigenvector=LHAPDF::xfx(xval,sqrt(qval),i);
		}
		else {
		  LHAPDF::initPDF( e+1 );
		  eigenvector=LHAPDF::xfx(xval,sqrt(qval),i);
		}
		val=val + (eigenvector - centralvalue_FS0) * abs(r);
		//  val=val + 0.5*abs(eigenvector1 - eigenvector) * (rand);
		// if (rand>0) { 
                //   LHAPDF::initPDF( e+1 );
		//   eigenvector=LHAPDF::xfx(xval,sqrt(qval),i);
		// }
		// else {
		//   LHAPDF::initPDF( e );
		//   eigenvector=LHAPDF::xfx(xval,sqrt(qval),i);
		// }
		// val=val + (eigenvector - centralvalue_FS0) * abs(rand);
              }        }  	     }
          else {
            for (int e=1; e<LHAPDF::numberPDF()+1; e++) { // symmetric eigenvector sets
	      double r = rand[(n-1)*LHAPDF::numberPDF()+e];
              LHAPDF::initPDF( e );
              double eigenvector=LHAPDF::xfx(xval,sqrt(qval),i);
              val=val+ (eigenvector - centralvalue_FS0) * (r); 
            }	  }
          if (abs(val)<1e-10) val=0;

          xfxval[0][x][q][(i+6)]=val;
          xfxval[1][x][q][(i+6)]=xfxval[1][x][q][(i+6)]+val;

          lhaout_errorsets << " "<<xfxval[0][x][q][(i+6)]<<"  "; 
        }
	      lhaout_errorsets <<endl; 
           
      }
    }
  }
 

  for (size_t x=0; x<npx; x++)
	{
	  for (size_t q=0; q<npq2; q++)
    {
      double xval=xg[x], qval=qg[q];
      for (int i=-6; i<7; i++) 
      {
        lhaout_centralval << " "<<xfxval[1][x][q][(i+6)]/nRandomReplicas<<"  ";  
      }              
		  lhaout_centralval <<endl; 
    }
	}

  cout << "Writing central value PDF: ";
  cout<<"\n\033[F\033[J";

  lhaout_errorsets << "'End:'"<<endl; 
  
  cout <<"LHAPDF Writeout successful" <<endl<<endl;
  delete [] xfxval;
  
  lhaout_header.close();
  lhaout_centralval.close();
  lhaout_errorsets.close();

  ifstream *infile;
  ofstream lhaout(ofilename.str().c_str());
  string line;
    
  string files[3] = {"header", "centralval", "errorsets" };
  for (int i=0; i<3; i++)
  {
    infile = new ifstream( (outdir+"/"+files[i]+".txt").c_str() );
    while ( getline ( *infile, line, '\n' ) )
    {
      lhaout << line << endl;
    }
    infile->close();
    remove((outdir+"/"+files[i]+".txt").c_str());
  }
  lhaout.close();

  return;
}
*/
