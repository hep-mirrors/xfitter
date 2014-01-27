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
 
  //set xmin xmax
  int nx = pdfgraphs.begin()->graph->GetN();
  if (opts.xmin == -1 && opts.xmax == -1)
    {
      opts.xmin = pdfgraphs.begin()->graph->GetX()[0];
      opts.xmax = pdfgraphs.begin()->graph->GetX()[nx - 1];
    }

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
  TCanvas *cnv = new TCanvas(cnvname, "", opts.resolution, opts.resolution);
  cnv->cd();
  if (opts.logx)
    cnv->SetLogx();
  cnv->SetLeftMargin(lmarg);
  cnv->SetRightMargin(rmarg);
  cnv->SetTopMargin(tmarg);

  TMultiGraph * mg = new TMultiGraph(((string)cnvname + "_multigraph").c_str(), "");
  TMultiGraph * mg_lines = new TMultiGraph(((string)cnvname + "_multigraph_lines").c_str(), "");
  double mx = 0;
  double mn = 0;
  for (vector <gstruct>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    {
      //Prepare graph line borders
      Double_t val_x[(*it).graph->GetN()], val_y[(*it).graph->GetN()], val_high_y[(*it).graph->GetN()], val_low_y[(*it).graph->GetN()]; 

      for (int i = 0; i < (*it).graph->GetN(); i++)
	{
	  Double_t val =  (*it).graph->GetY()[i];
	  Double_t errhigh =  (*it).graph->GetErrorYhigh(i);
	  Double_t errlow =  (*it).graph->GetErrorYlow(i);

	  val_x[i] = (*it).graph->GetX()[i];
	  val_high_y[i] = val + errhigh;
	  val_low_y[i] = val - errlow;
	  val_y[i] = val;
	}

      TGraph *centr = new TGraph((*it).graph->GetN(), val_x, val_y);
      TGraph *high = new TGraph((*it).graph->GetN(), val_x, val_high_y);
      TGraph *low = new TGraph((*it).graph->GetN(), val_x, val_low_y);

      for (int i = 0; i < (*it).graph->GetN(); i++)
	{
	  double xi = high->GetX()[i];
	  double yi_h = high->GetY()[i];
	  double yi_l = low->GetY()[i];
	  if (xi >= opts.xmin && xi <= opts.xmax)
	    {
	      mx = max(mx, yi_h);
	      mn = min(mn, yi_l);
	    }
	}

      //set border features      
      centr->SetLineColor((*it).graph->GetLineColor());
      centr->SetLineStyle(1);
      centr->SetLineWidth(2);
      high->SetLineColor((*it).graph->GetLineColor());
      high->SetLineStyle(1);
      high->SetLineWidth(2);
      low->SetLineColor((*it).graph->GetLineColor());
      low->SetLineStyle(1);
      low->SetLineWidth(2);

      //add graphs
      mg->Add((*it).graph);
      mg_lines->Add(centr);
      mg_lines->Add(high);
      mg_lines->Add(low);
    }

  //graphical settings
  mg->SetTitle(((string)" ; x  ; x" + pdflabels[ipdf] + "(x,Q^{2})").c_str());

  mg->Draw("ALE3"); //need to draw with A option to create axis

  if (mx != 0 || mn != 0)
    {
      float delta = mx - mn;
      mx = mx + delta * 0.8;
      mn = mn - delta * 0.1;
      mg->SetMaximum(mx);
      mg->SetMinimum(mn);
    }
  
  mg->GetXaxis()->Set(nx, opts.xmin, opts.xmax);
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
  mg_lines->Draw("L");

  //create legend
  TLegend * leg = new TLegend(lmarg+0.05, 1-tmarg-0.05-pdfgraphs.size()*0.05, lmarg+0.5, 1-tmarg-0.01);
  leg->SetTextSize(txtsize);
  leg->SetTextFont(62);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  //  leg->AddEntry((TObject*)0, ((string)"x" + pdflabels[ipdf] + " - Q^{2} = " + q2str + " GeV^{2}").c_str(), "");
  leg->AddEntry((TObject*)0, ((string)"Q^{2} = " + q2str + " GeV^{2}").c_str(), "");

  for (vector <gstruct>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    if (opts.dobands)
      leg->AddEntry((*it).graph, (*it).label.c_str(), "lf");
    else
      leg->AddEntry((*it).graph, (*it).label.c_str(), "l");

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
  TMultiGraph * mg_ratio_lines = new TMultiGraph(((string)cnvname + "_multigraph_ratio_lines").c_str(), "");

  //create legend
  TLegend * leg = new TLegend(lmarg+0.15, 1-tmarg-0.05-pdfgraphs.size()*0.05, lmarg+0.5, 1-tmarg-0.01);
  leg->SetTextFont(62);
  leg->SetTextSize(txtsize);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry((TObject*)0, ((string)"Q^{2} = " + q2str + " GeV^{2}").c_str(), "");

  //prepare ratio graph
  vector <gstruct>::iterator fit = pdfgraphs.begin();
  for (vector <gstruct>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    {
      TGraphAsymmErrors *r = (TGraphAsymmErrors*)(*it).graph->Clone();

      r->SetName(((string)(*it).graph->GetName() + "_ratio").c_str());

      Double_t val_x[(*it).graph->GetN()], val_y[(*it).graph->GetN()], val_high_y[(*it).graph->GetN()], val_low_y[(*it).graph->GetN()]; 

      for (int i = 0; i < (*it).graph->GetN(); i++)
	{
	  double val =  (*it).graph->GetY()[i];
	  double ref =  (*fit).graph->GetY()[i];
	  if (ref != 0)
	    r->SetPoint(i, (*it).graph->GetX()[i], val/ref);
	  else
	    r->SetPoint(i, (*it).graph->GetX()[i], 1);

	  if (opts.relerror)
	    r->SetPoint(i, (*it).graph->GetX()[i], 1);

	  if (opts.abserror)
	    r->SetPoint(i, (*it).graph->GetX()[i], 0);

	  double errhigh =  (*it).graph->GetErrorYhigh(i);
	  double errlow =  (*it).graph->GetErrorYlow(i);

	  double rathigh, ratlow;
	  rathigh = ( ref != 0)? (errhigh/ref) : 0;
	  ratlow = ( ref != 0)? (errlow/ref) : 0;
	  if (opts.relerror)
	    {
	      rathigh = ( val != 0 )? (errhigh/val) : 0;
	      ratlow = ( val != 0 )? (errlow/val) : 0;
	    }
	  if (opts.abserror)
	    {
	      rathigh = errhigh - val;
	      ratlow = errlow - val;
	    }
	  
	  r->SetPointError(i, 0, 0, ratlow, rathigh);
	}

      //prepare borders
      for (int i = 0; i < r->GetN(); i++)
	{
	  Double_t val =  r->GetY()[i];
	  Double_t errhigh =  r->GetErrorYhigh(i);
	  Double_t errlow =  r->GetErrorYlow(i);

	  val_x[i] = (*it).graph->GetX()[i];
	  val_y[i] = val;
	  val_high_y[i] = val+errhigh;
	  val_low_y[i] = val-errlow;
	}

      TGraph *r_centr = new TGraph((*it).graph->GetN(), val_x, val_y);
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
      r_centr->SetLineColor(opts.colors[colindx]);
      r_centr->SetLineStyle(1);
      r_centr->SetLineWidth(2);
      r_high->SetLineColor(opts.colors[colindx]);
      r_high->SetLineStyle(1);
      r_high->SetLineWidth(2);
      r_low->SetLineColor(opts.colors[colindx]);
      r_low->SetLineStyle(1);
      r_low->SetLineWidth(2);

      //add graphs
      mg_ratio->Add(r);
      mg_ratio_lines->Add(r_centr);
      mg_ratio_lines->Add(r_high);
      mg_ratio_lines->Add(r_low);

      //Add legend
      if (opts.dobands)
	leg->AddEntry((*it).graph, (*it).label.c_str(), "lf");
      else
	leg->AddEntry((*it).graph, (*it).label.c_str(), "l");
      colindx++;
    }

  //Make the TCanvas
  TCanvas *cnv = new TCanvas(((string)cnvname + "_ratio").c_str(), "pdf", opts.resolution, opts.resolution);
  cnv->cd();
  if (opts.logx)
    cnv->SetLogx();
  cnv->SetLeftMargin(lmarg);
  cnv->SetRightMargin(rmarg);
  cnv->SetTopMargin(tmarg);

  //graphical settings
  //  mg_ratio->SetTitle(((string)" ; x  ; x" + pdflabels[ipdf] + "(x,Q^{2})/x" + pdflabels[ipdf] + "(x,Q^{2})_{ref}").c_str());
  mg_ratio->Draw("ALE3"); //need to draw with A option to create axis

  mg_ratio->GetXaxis()->SetTitle(" x  ");
  mg_ratio->GetYaxis()->SetTitle(((string)" x" + pdflabels[ipdf] + "(x,Q^{2})/x" + pdflabels[ipdf] + "(x,Q^{2})_{ref}").c_str());
  if (opts.relerror)
    mg_ratio->GetYaxis()->SetTitle(((string)" #deltax" + pdflabels[ipdf] + "/#deltax" + pdflabels[ipdf] + "_{ref}").c_str());
  if (opts.abserror)
    mg_ratio->GetYaxis()->SetTitle(((string)" #deltax" + pdflabels[ipdf] + "").c_str());

  mg_ratio->SetMaximum(opts.rmax);
  mg_ratio->SetMinimum(opts.rmin);

  int nx = pdfgraphs.begin()->graph->GetN();
  mg_ratio->GetXaxis()->Set(nx, opts.xmin, opts.xmax);
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
  mg_ratio_lines->Draw("L");

  leg->Draw();

  DrawLogo("dc")->Draw();

  return cnv;
}

