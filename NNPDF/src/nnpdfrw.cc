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
#include "LHAPDF/FortranWrappers.h"

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
  
  TCanvas *dCpa = new TCanvas("paPlot", "P(#alpha)",12,38,699,499);		
  TGraph* dpalpha= new TGraph(NPOINTS,alpha,palph);
    
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
	
  dCpa->Print(("output/"+paout).c_str());

  //  exit(0);
  return;
}


/// callable convolution wrapper: create PDF replicas
extern "C" void create_randompdfreplicas_( const char* pdfset, const char* directory, int * nrep, int  * isSymmetric) {

  int nRandomReplicas = *nrep;
  stringstream ofilename;
  ofilename.str("");
  ofilename << "output/" << directory << "/" << pdfset << Form("_%dInputReplicas.LHgrid",nRandomReplicas);
  std::cout << "nnpdfrw isSymmetric  " << *isSymmetric<< std::endl;
  int hasSymmErrorSets = *isSymmetric;

  ifstream file(ofilename.str().c_str());
  if(file){
    std::cout<<"\033[1;31mOutput file " <<ofilename.str()<<" already exists. Skipping generation of random replicas!\033[0m"<<std::endl;
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

  size_t npx=100, npq2=50;
  typedef double xdim[npx][npq2][13];
  
  // Init x and q2 grid points
  vector<double> xg, qg;
  double qmin=LHAPDF::getQ2min(0);
  double qmax=LHAPDF::getQ2max(0);
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
    qg.push_back(qmin*pow(qmax/qmin,(((double) i ) -1)/(((double) npq2 ) -1)));
  
  
  double oas = LHAPDF::getOrderAlphaS();
  double opdf = LHAPDF::getOrderPDF();
  
  cout <<endl<<"Writing out replica LHAPDF grid: "<< pdfset << "_InputReplicas.LHgrid" <<endl;
  if (hasSymmErrorSets) cout <<endl<<"creating "<< nRandomReplicas << " random replicas PDF with symmetric errors" <<endl;
  else  cout <<endl<<"creating "<< nRandomReplicas << " random replicas from PDF with asymmetric errors" <<endl;
  cout <<"Using LHAPDF version: "<< LHAPDF::getVersion()<<endl;



  cout << "Creating random replica input LHgrid file" << endl;
  cout << ofilename.str() << endl;
  cout << "------------------------------------------" << endl;

  string dir = string(directory);
  mkdir( ("output/"+dir).c_str() ,0777);

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
  for (size_t i=0; i<(nRandomReplicas+1); i++)
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
  xdim *xfxval = new xdim[2];

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
