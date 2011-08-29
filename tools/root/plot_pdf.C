// #include "AtlasStyle.h"
// #include "AtlasStyle.C"
#include "Rtypes.h"

void plot_mc(string name="g", double q2=1.9, Int_t n=0)
{
//    SetAtlasStyle();
    TGraphErrors *p;

    TString nn;
    nn.Form("ave_%s_vs_x_for_Q2=%g",name.c_str(),q2);
    gDirectory->GetObject(nn,p);

//    p->Print();    
    
    Double_t ratsize = 0.0;

    TCanvas *c = new TCanvas("PDF","pdf",600,600);
    TPad *pad1 = new TPad("pad1","pad1",0.,ratsize,1.,1.);
    pad1->SetLogx();
    pad1->Draw();
    pad1->cd();
    p->GetXaxis()->Set(101,0.0001,1.);
    p->SetFillColor(kRed-2);
    p->SetFillStyle(3001);
    p->SetLineWidth(1);
    p->SetLineColor(kRed);
    p->GetYaxis()->SetTitle(name.c_str());

    p->GetYaxis()->SetTitleSize(0.06);
    p->GetYaxis()->SetTitleOffset(1.);
    p->GetXaxis()->Set(101,0.0001,1.);
    p->Draw("ALE3");

    Double_t av = 0;
    Double_t av2 = 0;
    
    for (Int_t i = 1; i<n+1 ; i++) {
        nn.Form("%s_vs_x_for_Q^{2}=%g_%i",name.c_str(),q2,i);
        TGraph *pp = NULL;        
        gDirectory->GetObject(nn,pp);        
        if (pp != NULL) {
            pp->SetLineColor(kGreen);
            pp->Draw("L");            
            
            av  += pp->GetY()[0];
            av2 += pp->GetY()[0]*pp->GetY()[0];
            cout << i << " "<<pp->GetY()[0] << endl;
            
            
        }        
    }
    
    av  /= n;
    av2 = sqrt(av2/n - av*av);
//    cout << n << " "<<av << " "<< av2<<endl;

    
    p->Draw("E3C");
    
    
}


void plot_pdf(string name="g", double q2=1.9, string dir="output/Graphs.root",
	  string name2=" ",string dir2=" ")
{
    SetAtlasStyle();
    TFile t(dir.c_str());
    TGraphAsymmErrors *p;
    TString nn;
    nn.Form("%s_vs_x_for_Q^{2}=%g",name.c_str(),q2);
    gDirectory->GetObject(nn,p);
//    gDirectory->ls();
    

    if (p == NULL) {
        cout << "Can not find graph "<< nn.Data()<<endl;
        return;        
    }
//    p->Print();
    
    
    TGraphAsymmErrors *p2 = NULL; 
    TGraphAsymmErrors *r2 = NULL;
  
    if (name2 != " ") {
        if (dir2 != " ") {
            TFile t(dir2.c_str());
        }
        nn.Form("%s_vs_x_for_Q^{2}=%g",name2.c_str(),q2);
        gDirectory->GetObject(nn,p2);
        p2->SetFillStyle(3001);
        p2->SetFillColor(kBlue);
        r2 = (TGraphAsymmErrors*)p2->Clone();
    }
    
    TGraphAsymmErrors *r = (TGraphAsymmErrors*)p->Clone();
  
    for(Int_t i=0; i<p->GetN(); i++) {
        r->SetPoint(i, p->GetX()[i], 1.);
        
        Double_t err =  p->GetErrorY(i);
        Double_t val =  p->GetY()[i];
        Double_t rat = ( val != 0 )? err/val : 0;
        r->SetPointError(i, 0, 0, rat, rat);
        
        if ( p2 != NULL) {
//        r2->SetPoint(i, p2->GetX()[i],  p2->GetY()[i]/ p->GetY()[i]);
            r2->SetPoint(i, p2->GetX()[i],  1.);
        


            Double_t err =  p2->GetErrorY(i);
            
            Double_t val =  p2->GetY()[i];

            Double_t rat = ( val != 0 )? err/val : 0;

            r2->SetPointError(i, 0, 0, rat, rat);

        
        }
    
    }

    Double_t ratsize = 0.3;

    TCanvas *c = new TCanvas("PDF","pdf",600,600);
    TPad *pad1 = new TPad("pad1","pad1",0.,ratsize,1.,1.);
    pad1->SetLogx();
    pad1->Draw();
    pad1->cd();
    p->GetXaxis()->Set(101,0.0001,1.);
    p->SetFillColor(kRed-10);
    p->SetFillStyle(1001);
    p->SetLineWidth(1);
    p->SetLineColor(kRed-2);
    p->GetYaxis()->SetTitle(name.c_str());
    
    p->GetYaxis()->SetTitleSize(0.06);
    p->GetYaxis()->SetTitleOffset(1.);
    p->GetXaxis()->Set(101,0.0001,1.);
    
    p->Draw("ACE3");

    if (p2 != NULL) {
        p2->SetLineColor(kBlue);      
        p2->Draw("CE3");
    }

    c->cd();
    TPad *pad2 = new TPad("pad1","pad1",0.,0.,1.,ratsize);
    pad1->SetBottomMargin(0.);
    pad2->SetTopMargin(0.);
    pad2->SetBottomMargin(0.3);
    pad2->SetLogx();
    pad2->Draw();
    pad2->cd();
    r->SetMaximum(1.35);
    r->SetMinimum(0.65);
    r->GetYaxis()->SetNdivisions(504);
    r->GetYaxis()->SetLabelSize(0.12);
    r->GetXaxis()->SetLabelSize(0.12);
    r->GetXaxis()->SetLabelOffset(0.03);
    r->GetXaxis()->Set(101,0.0001,1.);
    r->GetXaxis()->SetTitle("x  ");
    r->GetXaxis()->SetTitleSize(0.12);
    r->GetXaxis()->SetTitleOffset(1.);
    
    //r->GetXaxis()->SetTitleOffset(0.5);  
    //  r->GetXaxis()->SetTitleSize(.15);  
    r->SetFillColor(2);
    r->SetFillStyle(1001);
    
    r->SetFillColor(kRed-10);
    r->SetLineStyle(1);
    r->SetLineWidth(1);
  
  
    r->Draw("ALE3");

    if (r2 != NULL) {
        r2->SetFillColor(kBlue);
        r2->SetFillStyle(3001);
        r2->Draw("E3");
    }
  
    
}

