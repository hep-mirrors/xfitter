#include <PdfsPainter.h>

#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TMath.h>

#include <iostream>

TCanvas * PdfsPainter(double q2, string pdf, vector <gstruct> pdfgraphs)
{
  if (pdfgraphs.size() < 1)
    {
      cout << "Empty pdf TGraph vector for pdf: " << pdf << endl;
      return 0;
    }

  char q2str[30];				
  if (q2 < 10)
    sprintf(q2str, "%.1f",  q2);
  else
    sprintf(q2str, "%.0f",  q2);
  
  char cnvname[30];
  sprintf(cnvname, "q2_%s_pdf_%s",  q2str, pdf.c_str());


  //prepare TGraphs
  int colors[] = {kRed, kBlue, kMagenta, (kGreen+2), kCyan, kYellow};
  int styles[] = {3004, 3005, 3006, 3007, 3016, 3020};
  int colindx = 0;
 
  //set graph features
  for (vector <gstruct>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    {
      (*it).graph->GetXaxis()->Set(101,0.0001,1.);
      (*it).graph->SetFillColor(colors[colindx]);
      (*it).graph->SetFillStyle(styles[colindx]);
      (*it).graph->SetLineStyle(1);
      (*it).graph->SetLineWidth(1);
      (*it).graph->SetLineColor(colors[colindx]);
      (*it).graph->GetYaxis()->SetTitleSize(0.06);
      (*it).graph->GetYaxis()->SetTitleOffset(1.);

      colindx++;
    }	  

  //Make the TCanvas
  TCanvas *cnv = new TCanvas(cnvname, "", 800, 600);
  cnv->cd();
  cnv->SetLogx();

  TMultiGraph * mg = new TMultiGraph(((string)cnvname + "_multigraph").c_str(), "");
  for (vector <gstruct>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    {
      (*it).graph->SetMaximum(TMath::MaxElement((*it).graph->GetN(), (*it).graph->GetY()));
      float mx = (*it).graph->GetMaximum();
      (*it).graph->SetMinimum(TMath::MinElement((*it).graph->GetN(), (*it).graph->GetY()));
      float mn = (*it).graph->GetMinimum();     
      mx = mx + (mx - mn) * 0.8;
      //     mn = mn - (mx - mn) * 0.3;       
      mg->SetMaximum(mx);
      //     mg->SetMinimum(mn);

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
      r_high->SetLineWidth(1);
      r_low->SetLineColor((*it).graph->GetLineColor());
      r_low->SetLineStyle(1);
      r_low->SetLineWidth(1);

      //add graphs
      mg->Add((*it).graph);
      mg->Add(r_high);
      mg->Add(r_low);
    }

  //graphical settings
  mg->SetTitle(" ; x  ; xf(x,Q^{2})");

  mg->Draw("ALE3"); //need to draw with A option to create axis

  mg->GetXaxis()->Set(101,0.0001,1.);    
  mg->GetXaxis()->SetTitleFont(62);
  mg->GetXaxis()->SetLabelFont(62);
  mg->GetXaxis()->SetTitleSize(0.05);
  mg->GetXaxis()->SetLabelSize(0.05);
  mg->GetXaxis()->SetTitleOffset(1.);

  mg->GetYaxis()->SetTitleFont(62);
  mg->GetYaxis()->SetLabelFont(62);
  mg->GetYaxis()->SetTitleSize(0.05);
  mg->GetYaxis()->SetLabelSize(0.05);      
  mg->GetYaxis()->SetTitleOffset(1.);

  mg->Draw("ALE3");

  //create legend
  TLegend * leg = new TLegend(0.15, 0.7, 0.45, 0.9);
  leg->AddEntry((TObject*)0, ((string)"x" + pdf + "  at  " + "Q^{2} = " + q2str + " GeV^{2}").c_str(), "");

  for (vector <gstruct>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    {
      leg->AddEntry((*it).graph, (*it).label.c_str());
    }
  leg->SetTextSize(0.048);
  leg->SetTextFont(62);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->Draw();

  return cnv;
}

TCanvas * PdfsRatioPainter(double q2, string pdf, vector <gstruct> pdfgraphs)
{
    if (pdfgraphs.size() < 1)
    {
      cout << "Empty pdf TGraph vector for pdf: " << pdf << endl;
      return 0;
    }

  char q2str[30];				
  if (q2 < 10)
    sprintf(q2str, "%.1f",  q2);
  else
    sprintf(q2str, "%.0f",  q2);
  
  char cnvname[30];
  sprintf(cnvname, "q2_%s_pdf_%s",  q2str, pdf.c_str());

  //prepare TGraphs
  int colors[] = {kRed, kBlue, kMagenta, (kGreen+2), kCyan, kYellow};
  int styles[] = {3004, 3005, 3006, 3007, 3016, 3020};
  int colindx = 0;
  TMultiGraph * mg_ratio = new TMultiGraph(((string)cnvname + "_multigraph_ratio").c_str(), "");

  //create legend
  TLegend * leg = new TLegend(0.15, 0.7, 0.45, 0.9);
  leg->AddEntry((TObject*)0, (pdf + "  at  Q^{2} = " + q2str + " GeV^{2}").c_str(), "");

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

      //Graphical settings
      r->SetMaximum(1.35);
      r->SetMinimum(0.65);

      r->GetXaxis()->SetTitle("x  ");
      r->GetXaxis()->SetTitleSize(0.12);
      r->GetXaxis()->SetTitleOffset(1.);
      r->GetXaxis()->Set(101,0.0001,1.);
      r->GetXaxis()->SetLabelOffset(0.03);
      r->GetXaxis()->SetLabelSize(0.12);

      r->GetYaxis()->SetTitle("pdf(g,u,d,..)");
      r->GetYaxis()->SetTitleSize(0.06);
      r->GetYaxis()->SetTitleOffset(1.);
      r->GetYaxis()->SetNdivisions(504);
      r->GetYaxis()->SetLabelSize(0.12);
      
      //set graph features
      r->SetFillColor(colors[colindx]);
      r->SetFillStyle(styles[colindx]);
      r->SetLineColor(colors[colindx]);
      r->SetLineStyle(1);
      r->SetLineWidth(1);

      //set border features
      r_high->SetLineColor(colors[colindx]);
      r_high->SetLineStyle(1);
      r_high->SetLineWidth(1);
      r_low->SetLineColor(colors[colindx]);
      r_low->SetLineStyle(1);
      r_low->SetLineWidth(1);

      //add graphs
      mg_ratio->Add(r);
      mg_ratio->Add(r_high);
      mg_ratio->Add(r_low);

      //set legend features
      leg->AddEntry((*it).graph, (*it).label.c_str());
      leg->SetTextFont(62);
      leg->SetTextSize(0.035);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);

      colindx++;
    }	  

  //Make the TCanvas
  TCanvas *cnv = new TCanvas(((string)cnvname + "_ratio").c_str(), "pdf", 800, 800);
  cnv->cd();
  cnv->SetLogx();

  //graphical settings
  mg_ratio->SetTitle(" ; x  ; xf(x,Q^{2})/xf(x,Q^{2})_{ref}");

  mg_ratio->Draw("ALE3"); //need to draw with A option to create axis

  mg_ratio->SetMaximum(3);
  mg_ratio->SetMinimum(-1);

  mg_ratio->GetXaxis()->SetTitleFont(62);
  mg_ratio->GetXaxis()->SetLabelFont(62);
  mg_ratio->GetXaxis()->SetTitleSize(0.035);
  mg_ratio->GetXaxis()->SetLabelSize(0.035);
  mg_ratio->GetXaxis()->SetTitleOffset(1.4);
  
  mg_ratio->GetYaxis()->SetTitleFont(62);
  mg_ratio->GetYaxis()->SetLabelFont(62);
  mg_ratio->GetYaxis()->SetTitleSize(0.035);
  mg_ratio->GetYaxis()->SetLabelSize(0.035);
  mg_ratio->GetYaxis()->SetTitleOffset(1.4);

  leg->Draw();

  return cnv;
}

