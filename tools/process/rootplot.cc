#include <TH1F.h>
#include <TCanvas.h>
#include <stdlib.h>
#include "utils.h"
#include <string.h>
#include <libgen.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <TAxis.h>
#include <TGraph.h>
#include <TStyle.h>
#include <numeric>

using namespace std;

extern "C" void testplot();
//extern "C" double Palpha(double * chi2_alpha, double alpha, int nchi2);
extern "C" void FillwHist(double * weights,  double * chi2,int ndata,double nrep_old, int method);
extern "C" double calc_Palpha(double alpha) ;

void testplot(){
  return;
};

void BinLogX(TH1*h)
{
        
    TAxis *axis = h->GetXaxis();
    int bins = axis->GetNbins();
        
    Axis_t from = axis->GetXmin();
    Axis_t to = axis->GetXmax();
    Axis_t width = (to - from) / bins;
    Axis_t *new_bins = new Axis_t[bins + 1];
        
    for (int i = 0; i <= bins; i++) {
        new_bins[i] = pow(10, from + i * width);
    }
    axis->Set(bins, new_bins);
    delete new_bins;
}

void FillwHist(double * weights,  double * chi2,int ndata, double nrep_old, int method){

	// Weight Histogram
	TCanvas *wHc = new TCanvas ("wHistC","Weight Histogram",12,38,699,499);
	TH1F* wHist= new TH1F("wH", "Weight Histogram;Weight;Frequency;", 50, -7, 1);
  
	// Log binning  
	BinLogX(wHist);
	wHc->GetPad(0)->SetLogx();

       for(int i=1; i< nrep_old; i++) { 
	 if (weights[i]!=0) {
	   wHist->Fill(weights[i]);
	 }
       }

	//// TO DO CHECK PLOTS
	 wHist->Draw(); 
	 wHc->Print(("./weights.pdf")); 
	 delete wHc; 

	 vector<double> logw;  
	 vector<double> w;  

	 // P(alpha) plot
  	     double p=0;
 	     double ptot=0;

	 size_t NPOINTS=100;
	 Double_t alpha[NPOINTS],palph[NPOINTS];
	 // double ptot=0;
	 for (size_t np=0; np<NPOINTS; np++)
	   {
	     logw.clear();
	     w.clear();

	     alpha[np]=5*(Double_t ) np/(Double_t  )NPOINTS + 0.1; 

	     double w_alpha[(int)nrep_old];
	     double chi2_alpha[(int)nrep_old];  //chi2;
  
	     /// calculate the palpha here
	     for (size_t j=0; j<nrep_old; j++) {
	       chi2_alpha[j]=(chi2[j])/(alpha[np]*alpha[np]);
	     }
	       
	     /// now compute weights
	     for(int i=0; i< nrep_old; i++) { 
	       if (method==0) logw.push_back( -1* (chi2_alpha[i-1])/(2.0)); //Giele-Keller
	       else   logw.push_back( - chi2_alpha[i-1]/(2.0) +( (((double) ndata)-1.0)/2.)*log(chi2_alpha[i-1])); //Bayesian
	     }

	     // Calculate weights
	     for (size_t i=0 ;i<nrep_old; i++){
	       w.push_back(exp(logw[i] ));
	     }

	     double wtot = accumulate(w.begin(),w.end(),0.0); 
	     for (size_t i=0; i<w.size(); i++) {
	       if ((w[i]*(nrep_old/wtot)) < 1e-12)
		 w[i]=0;
	     }

	     p=0;
	     for (size_t j=0; j<w.size(); j++) {
	       p+=w[j];
	     }
	     
	     palph[np]=p/alpha[np];
	          
	     ptot=ptot+palph[np];
	         
	   }
	
	 // Roughly normalise
	 double integral=0;
    
	 integral+=palph[0]+palph[NPOINTS-1];
    
	 for ( size_t j=1; j<(NPOINTS)/2 ; j++ )
	   integral+=2*palph[2*j -1];
 
	 for (size_t j=1; j<(NPOINTS)/2 + 1; j++)
	   integral+=4*palph[2*j - 2];
	 
	 double h = 5.0/(double) NPOINTS;
	 integral= integral*h/3.0;

	 
	 for (size_t i=0; i<NPOINTS; i++)
	   palph[i]=palph[i]/(integral);

	 TCanvas *dCpa = new TCanvas("paPlot", "P(alpha)",12,38,699,499);	
	 TGraph* dpalpha= new TGraph(NPOINTS,alpha,palph);
	 double alphamin=5*(Double_t ) 0/(Double_t  )NPOINTS + 0.1; 
	 double alphamax= 5*(Double_t ) (NPOINTS-1)/(Double_t  )NPOINTS + 0.1; 
	 TH1F* hCpa= new TH1F("hCpa", "pa Plot;(#alpha);P(#alpha);", 100,alphamin , alphamax);
 
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
    
	 hCpa->SetMaximum(2.2);
	 hCpa->Draw();
	 dpalpha->Draw("L");
	 if (method==1) 	dCpa->Print(("./palpha.pdf")); 
	 
	 return;
}

