#include <PdfsPainter.h>
#include <DrawLogo.h>
#include <CommandParser.h>

#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TMath.h>

#include <iostream>

//useful trick to declare a vector "statically"
//be carefull! pdf names must be mapped to pdf definition in Output.h
//  enum pdf{kNULL=-1, kGluon=0, kU=1, kD=2, kUv=3, kDv=4, kUb=5, kDb=6, kSea=7, kUSea=8, kDSea=9, kS=10, kC=11, kB=12};
string pdflab[] = {"g", "U", "D", "uv", "dv", "#bar{U}", "#bar{D}", "sea", "#bar{u}", "#bar{d}", "s", "c", "b"};
string pdffil[] = {"g", "U", "D", "uv", "dv", "UBar", "DBar", "sea", "ubar", "dbar", "s", "c", "b"};

vector <string> pdflabels(pdflab, pdflab + sizeof(pdflab) / sizeof(string));
vector <string> pdffiles(pdffil, pdffil + sizeof(pdffil) / sizeof(string));


float txtsize = 0.045;
float offset = 1.5;
float lmarg = 0.15;
//float bmarg = 0.2;

TCanvas * PdfsPainter(double q2, int ipdf, vector <gstruct> pdfgraphs)
{
  if (pdfgraphs.size() < 1)
    {
      cout << "Empty pdf TGraph vector for pdf: " << pdffiles[ipdf] << endl;
      return 0;
    }

  char q2str[30];				
  if (q2 < 10)
    sprintf(q2str, "%.1f",  q2);
  else
    sprintf(q2str, "%.0f",  q2);
  
  char cnvname[30];
  sprintf(cnvname, "q2_%s_pdf_%s",  q2str, pdffiles[ipdf].c_str());


  //prepare TGraphs
  int colindx = 0;
 
  //set graph features
  for (vector <gstruct>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    {
      (*it).graph->SetFillColor(opts.colors[colindx]);
      if (opts.filledbands)
	(*it).graph->SetFillStyle(1001);
      else
	(*it).graph->SetFillStyle(opts.styles[colindx]);
      (*it).graph->SetLineStyle(1);
      (*it).graph->SetLineWidth(2);
      (*it).graph->SetLineColor(opts.colors[colindx]);
      colindx++;
    }	  

  //Make the TCanvas
  TCanvas *cnv = new TCanvas(cnvname, "", 2400, 2400);
  cnv->cd();
  cnv->SetLogx();
  cnv->SetLeftMargin(lmarg);

  TMultiGraph * mg = new TMultiGraph(((string)cnvname + "_multigraph").c_str(), "");
  double mx = 0;
  double mn = 0;
  for (vector <gstruct>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    {
      (*it).graph->SetMaximum(TMath::MaxElement((*it).graph->GetN(), (*it).graph->GetY()));
      mx = max(mx,(*it).graph->GetMaximum());
      (*it).graph->SetMinimum(TMath::MinElement((*it).graph->GetN(), (*it).graph->GetY()));
      mn = min(mn, (*it).graph->GetMinimum());

      //prepare borders
      Double_t val_x[(*it).graph->GetN()], val_high_y[(*it).graph->GetN()], val_low_y[(*it).graph->GetN()]; 

      for (int i = 0; i < (*it).graph->GetN(); i++)
	{
	  Double_t val =  (*it).graph->GetY()[i];
	  Double_t errhigh =  (*it).graph->GetErrorYhigh(i);
	  Double_t errlow =  (*it).graph->GetErrorYlow(i);

	  val_x[i] = (*it).graph->GetX()[i];
	  val_high_y[i] = val+errhigh;
	  val_low_y[i] = val-errlow;
	}

      TGraph *r_high = new TGraph((*it).graph->GetN(), val_x, val_high_y);
      TGraph *r_low = new TGraph((*it).graph->GetN(), val_x, val_low_y);

      //set border features      
      r_high->SetLineColor((*it).graph->GetLineColor());
      r_high->SetLineStyle(1);
      r_high->SetLineWidth(2);
      r_low->SetLineColor((*it).graph->GetLineColor());
      r_low->SetLineStyle(1);
      r_low->SetLineWidth(2);

      //add graphs
      mg->Add((*it).graph);
      mg->Add(r_high);
      mg->Add(r_low);
    }

  //graphical settings
  mg->SetTitle(" ; x  ; xf(x,Q^{2})");

  mg->Draw("ALE3"); //need to draw with A option to create axis

  if (mx != 0 || mn != 0)
    {
      mx = mx + (mx - mn) * 0.8;
      //     mn = mn - (mx - mn) * 0.3;       
      mg->SetMaximum(mx);
      //     mg->SetMinimum(mn);
    }

  mg->GetXaxis()->Set(101,0.0001,1.);    
  mg->GetXaxis()->SetTitleFont(62);
  mg->GetXaxis()->SetLabelFont(62);
  mg->GetXaxis()->SetTitleSize(txtsize);
  mg->GetXaxis()->SetLabelSize(txtsize);
  //  mg->GetXaxis()->SetTitleOffset(offset);

  mg->GetYaxis()->SetTitleFont(62);
  mg->GetYaxis()->SetLabelFont(62);
  mg->GetYaxis()->SetTitleSize(txtsize);
  mg->GetYaxis()->SetLabelSize(txtsize);      
  mg->GetYaxis()->SetTitleOffset(offset);

  mg->Draw("ALE3");

  //create legend
  TLegend * leg = new TLegend(0.15, 0.7, 0.45, 0.9);
  leg->AddEntry((TObject*)0, ((string)"x" + pdflabels[ipdf] + "  at  " + "Q^{2} = " + q2str + " GeV^{2}").c_str(), "");

  for (vector <gstruct>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    {
      leg->AddEntry((*it).graph, (*it).label.c_str());
    }
  leg->SetTextSize(txtsize);
  leg->SetTextFont(62);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->Draw();

  DrawLogo()->Draw();

  return cnv;
}

TCanvas * PdfsRatioPainter(double q2, int ipdf, vector <gstruct> pdfgraphs)
{
    if (pdfgraphs.size() < 1)
    {
      cout << "Empty pdf TGraph vector for pdf: " << pdffiles[ipdf] << endl;
      return 0;
    }

  char q2str[30];				
  if (q2 < 10)
    sprintf(q2str, "%.1f",  q2);
  else
    sprintf(q2str, "%.0f",  q2);
  
  char cnvname[30];
  sprintf(cnvname, "q2_%s_pdf_%s",  q2str, pdffiles[ipdf].c_str());

  //prepare TGraphs
  int colindx = 0;
  TMultiGraph * mg_ratio = new TMultiGraph(((string)cnvname + "_multigraph_ratio").c_str(), "");

  //create legend
  TLegend * leg = new TLegend(0.15, 0.7, 0.45, 0.9);
  leg->AddEntry((TObject*)0, (pdflabels[ipdf] + "  at  Q^{2} = " + q2str + " GeV^{2}").c_str(), "");

  //prepare ratio graph
  vector <gstruct>::iterator fit = pdfgraphs.begin();
  for (vector <gstruct>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    {
      TGraphAsymmErrors *r = (TGraphAsymmErrors*)(*it).graph->Clone();

      r->SetName(((string)(*it).graph->GetName() + "_ratio").c_str());

      Double_t val_x[(*it).graph->GetN()], val_high_y[(*it).graph->GetN()], val_low_y[(*it).graph->GetN()]; 

      for (int i = 0; i < (*it).graph->GetN(); i++)
	{
	  Double_t val =  (*it).graph->GetY()[i];
	  Double_t ref =  (*fit).graph->GetY()[i];
	  r->SetPoint(i, (*it).graph->GetX()[i], val/ref);

	  Double_t errhigh =  (*it).graph->GetErrorYhigh(i);
	  Double_t errlow =  (*it).graph->GetErrorYlow(i);
	  Double_t rathigh = ( val != 0 )? ((val+errhigh)/ref-val/ref) : 0;
	  Double_t ratlow = ( val != 0 )? ((val+errlow)/ref-val/ref) : 0;
	  r->SetPointError(i, 0, 0, ratlow, rathigh);
	}

      //prepare borders
      for (int i = 0; i < r->GetN(); i++)
	{
	  Double_t val =  r->GetY()[i];
	  Double_t errhigh =  r->GetErrorYhigh(i);
	  Double_t errlow =  r->GetErrorYlow(i);

	  val_x[i] = (*it).graph->GetX()[i];
	  val_high_y[i] = val+errhigh;
	  val_low_y[i] = val-errlow;
	}

      TGraph *r_high = new TGraph((*it).graph->GetN(), val_x, val_high_y);
      TGraph *r_low = new TGraph((*it).graph->GetN(), val_x, val_low_y);

      //set graph features
      r->SetFillColor(opts.colors[colindx]);
      if (opts.filledbands)
	r->SetFillStyle(1001);
      else
	r->SetFillStyle(opts.styles[colindx]);
      r->SetLineColor(opts.colors[colindx]);
      r->SetLineStyle(1);
      r->SetLineWidth(2);

      //set border features
      r_high->SetLineColor(opts.colors[colindx]);
      r_high->SetLineStyle(1);
      r_high->SetLineWidth(2);
      r_low->SetLineColor(opts.colors[colindx]);
      r_low->SetLineStyle(1);
      r_low->SetLineWidth(2);

      //add graphs
      mg_ratio->Add(r);
      mg_ratio->Add(r_high);
      mg_ratio->Add(r_low);

      //set legend features
      leg->AddEntry((*it).graph, (*it).label.c_str());
      leg->SetTextFont(62);
      leg->SetTextSize(txtsize);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);

      colindx++;
    }	  

  //Make the TCanvas
  TCanvas *cnv = new TCanvas(((string)cnvname + "_ratio").c_str(), "pdf", 2400, 2400);
  cnv->cd();
  cnv->SetLogx();
  cnv->SetLeftMargin(lmarg);

  //graphical settings
  mg_ratio->SetTitle(" ; x  ; xf(x,Q^{2})/xf(x,Q^{2})_{ref}");

  mg_ratio->Draw("ALE3"); //need to draw with A option to create axis

  mg_ratio->SetMaximum(opts.rmax);
  mg_ratio->SetMinimum(opts.rmin);

  mg_ratio->GetXaxis()->Set(101,0.0001,1.);    
  mg_ratio->GetXaxis()->SetTitleFont(62);
  mg_ratio->GetXaxis()->SetLabelFont(62);
  mg_ratio->GetXaxis()->SetTitleSize(txtsize);
  mg_ratio->GetXaxis()->SetLabelSize(txtsize);
  //  mg_ratio->GetXaxis()->SetTitleOffset(offset);
  
  mg_ratio->GetYaxis()->SetTitleFont(62);
  mg_ratio->GetYaxis()->SetLabelFont(62);
  mg_ratio->GetYaxis()->SetTitleSize(txtsize);
  mg_ratio->GetYaxis()->SetLabelSize(txtsize);
  mg_ratio->GetYaxis()->SetTitleOffset(offset);
  mg_ratio->GetYaxis()->SetNdivisions(504);

  mg_ratio->Draw("ALE3");

  leg->Draw();

  //  DrawLogo()->SetPad(0.64, 0.75, 0.79, 0.89);
  DrawLogo()->Draw();

  return cnv;
}

