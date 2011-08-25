// #include "AtlasStyle.h"
// #include "AtlasStyle.C"

void plot_pdf(string name="g", double q2=1.9, string dir="output/Graphs.root",
	  string name2=" ",string dir2=" "){
  //  SetAtlasStyle();
  TFile t(dir.c_str());
  TGraphAsymmErrors *p;
  TString nn;
  nn.Form("%s_vs_x_for_Q^{2}=%g",name.c_str(),q2);
  gDirectory->GetObject(nn,p);

  TGraphAsymmErrors *p2 = NULL; 
  if (name2 != " ") {
    if (dir2 != " ") {
      TFile t(dir2.c_str());
    }
    nn.Form("%s_vs_x_for_Q^{2}=%g",name2.c_str(),q2);
    gDirectory->GetObject(nn,p2);
    p2->SetFillStyle(3003);
    p2->SetFillColor(4);
  }
  
  TGraphAsymmErrors *r = (TGraphAsymmErrors*)p->Clone();

  cout << "AH"<<endl;
  for(Int_t i=0; i<p->GetN(); i++) {
    r->SetPoint(i, p->GetX()[i], 1.);

    Double_t err =  p->GetErrorY(i);
    Double_t val =  p->GetY()[i];
    Double_t rat = ( val != 0 )? err/val : 0;
    r->SetPointError(i, 0, 0, rat, rat);
  }

  Double_t ratsize = 0.3;

  TCanvas *c = new TCanvas("PDF","pdf",600,600);
  TPad *pad1 = new TPad("pad1","pad1",0.,ratsize,1.,1.);
  pad1->SetLogx();
  pad1->Draw();
  pad1->cd();
  p->GetXaxis()->Set(101,0.0001,1.);
  p->SetFillColor(2);
  p->SetFillStyle(3001);
  p->Draw("ACE3");

  if (p2 != NULL) {
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
  //r->GetXaxis()->SetTitleOffset(0.5);  
  //  r->GetXaxis()->SetTitleSize(.15);  
  r->SetFillColor(2);
  r->SetFillStyle(3001);
  r->Draw("ALE3"); 
}

