#include <stdlib.h>
#include <iostream>
#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "PdfTable.h"
#include <getopt.h>
#include "TAxis.h"
#include <math.h>
using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char **argv) {

    TFile output("Graphs.root","RECREATE");
    
    // Loop over all input files
    for (Int_t iq2=0; iq2<20; iq2++) {

        vector< vector<double> > tab1;
        tab1.clear();

        vector< vector<double> > tab2;
        tab2.clear();

        int npdf = 0;
        int nx   = 0;

        vector<string> names;
        names.clear();

        vector<double> xvals;
        xvals.clear();
        
        Double_t Q2 = 0;        

        
        for (int ifile=1; ifile<argc; ifile++) {
            TString filename("");
            filename.Form("%s/pdfs_q2val_%02d.txt",argv[ifile], iq2+1);
            PdfTable* table =  new PdfTable(filename.Data());

            if (table->GetNx() == 0)  break;

            nx = table->GetNx();
            
            
            npdf = table->GetNPdfs();
            Q2 = table->GetQ2();            
            
            if (ifile == 1) {
                cout << "Q2="<<Q2<<" npdf="<<npdf<<endl;
            
                xvals = table->GetPDF("x");
                // Prepare graphs with mean and RMS values
                for (Int_t i=1; i<=npdf; i++) {
                    string name = table->GetColumnName(i);

                    names.push_back(name);                    

                    vector<double> pdf = table->GetPDF(name);
                    tab1.push_back(pdf);                    
                    // square:
                    for (Int_t i =0; i< pdf.size(); i++) {
                        pdf[i] *= pdf[i];                        
                    }                
                    tab2.push_back(pdf);                    
                }
            }
            
                
            
            // Store all individual replica
            for (Int_t i=1; i<=npdf; i++) {

                string name = table->GetColumnName(i);
                TGraph *graph = table->GetPDFGraph(name);
                TString GraphName;
                GraphName.Form("%s_vs_x_for_Q^{2}=%g_%i",name.c_str(),Q2,ifile);
                TString GraphTitle;
                GraphTitle.Form("%s vs x for Q^{2}=%g GeV^{2}",name.c_str(),Q2);
                
                graph->SetName(GraphName);
                graph->SetTitle(GraphTitle);
                graph->GetXaxis()->SetTitle("x");
                graph->GetYaxis()->SetTitle(name.c_str());
                
                graph->Write();


                // Mean/RMS:
                if (ifile>1) {
                    
                    vector<double> pdf = table->GetPDF(name);
                
                    for (Int_t ii =0; ii< pdf.size(); ii++) {
//                    cout << name <<  " " << filename.Data()<< " "<< pdf.size() << endl;                    
//                    cout << "i"<<i<< " " << ii << " "<< tab1[i-1][ii]<<endl;
                    
                        tab1[i-1][ii] += pdf[ii];                        
                        tab2[i-1][ii] += pdf[ii]*pdf[ii];                        

//                        if ( (i-1 == 0 ) && (ii==0) ) {
//                            cout << pdf[ii];                            
//                        }
                        
                    }                
                    pdf.clear();
                }
                
                
            }            
            delete table;
        }
        // Get averages:
        for (Int_t i=0; i<npdf; i++) {
            
            TGraphErrors *Ave = new TGraphErrors(nx);
            for (Int_t j = 0; j<nx; j++) {
                tab1[i][j] /= (argc-1);
                tab2[i][j] = sqrt(fabs(tab2[i][j]/(argc-1)-tab1[i][j]*tab1[i][j]));                
                Ave->SetPoint(j,xvals[j],tab1[i][j]);
                Ave->SetPointError(j,0.,tab2[i][j]);
                
            }

            TString GraphName;
            GraphName.Form("ave_%s_vs_x_for_Q2=%g",names[i].c_str(),Q2);
            TString GraphTitle;
            GraphTitle.Form("ave %s vs x for Q^{2}=%g GeV^{2}",names[i].c_str(),Q2);

            Ave->SetName(GraphName.Data());
            Ave->SetTitle(GraphTitle.Data());

            Ave->Write();
            
            
        }
        
    }

    output.Write();
    output.Close();

        
    return 0;
}
