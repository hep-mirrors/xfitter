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

#include "Outdir.h"

vector <TCanvas*> PdfsPainter(double q2, pdftype ipdf)
{
  vector <TCanvas*> cnvs;

  vector <TGraphAsymmErrors*> pdfgraphs;
  vector <string> labels;
  for (vector<string>::iterator itl = opts.labels.begin(); itl != opts.labels.end(); itl++)
    if (pdfmap[*itl].Central.find(q2) !=  pdfmap[*itl].Central.end())
      {
	pdfgraphs.push_back(pdfmap[*itl].Central[q2].GetPdf(ipdf));
	labels.push_back(*itl);
      }

  if (pdfgraphs.size() < 1)
    {
      cout << "Error: Empty pdf TGraph vector for pdf: " << pdffiles[(int)ipdf] << endl;
      exit(1);
    }

  char q2str[30];				
  if (q2 < 10)
    sprintf(q2str, "%.1f",  q2);
  else
    sprintf(q2str, "%.0f",  q2);
  
  char cnvname[30];
  sprintf(cnvname, "q2_%s_pdf_%s",  q2str, pdffiles[ipdf].c_str());

  //Set xmin xmax
  for (vector <TGraphAsymmErrors*>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    {
      if (*it == 0)
	continue;
      if (opts.xmin == -1 && opts.xmax == -1)
	{
	  opts.xmin = (*it)->GetX()[0];
	  opts.xmax = (*it)->GetX()[(*it)->GetN() - 1];
	}
    }

  //Remove points out of x range
  for (vector <TGraphAsymmErrors*>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    {
      if (*it == 0)
	continue;
      bool removed = true;
      while (removed)
	{
	  removed = false;
	  for (int i = 0; i < (*it)->GetN(); i++)
	    {
	      double xi = (*it)->GetX()[i];
	      if (xi <= opts.xmin || xi >= opts.xmax)
		{
		  (*it)->RemovePoint(i);
		  removed = true;
		  break;
		}
	    }
	}
    }

  //Set colors and styles
  int colindx = 0;
  for (vector <TGraphAsymmErrors*>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    {
      if (*it == 0)
	continue;
      (*it)->SetFillColor(opts.colors[labels[it-pdfgraphs.begin()]]);
      if (opts.filledbands)
	(*it)->SetFillStyle(1001);
      else
	(*it)->SetFillStyle(opts.styles[labels[it-pdfgraphs.begin()]]);
      (*it)->SetLineStyle(1);
      (*it)->SetLineWidth(opts.lwidth);
      (*it)->SetLineColor(opts.colors[labels[it-pdfgraphs.begin()]]);
      colindx++;
    }	  

  //Prepare TGraphs
  TMultiGraph * mg = new TMultiGraph(((string)cnvname + "_multigraph").c_str(), "");
  TMultiGraph * mg_lines = new TMultiGraph(((string)cnvname + "_multigraph_lines").c_str(), "");
  TMultiGraph * mg_shade = new TMultiGraph(((string)cnvname + "_multigraph_shade").c_str(), "");
  double mx = 0;
  double mn = 0;
  for (vector <TGraphAsymmErrors*>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    {
      if (*it == 0)
	continue;
      //Prepare graph line borders and graph shade
      int npoints = (*it)->GetN();
      double val_x[npoints], val_y[npoints], val_high_y[npoints], val_low_y[npoints]; 
      double xsh[2*npoints], ysh[2*npoints];

      for (int i = 0; i < (*it)->GetN(); i++)
	{
	  double val = (*it)->GetY()[i];
	  double errhigh = (*it)->GetErrorYhigh(i);
	  double errlow = (*it)->GetErrorYlow(i);

	  val_x[i] = (*it)->GetX()[i];
	  val_y[i] = val;
	  val_high_y[i] = val + errhigh;
	  val_low_y[i] = val - errlow;

	  //shade TGraph
	  xsh[i] = (*it)->GetX()[i];
	  ysh[i] = val + errhigh;
	  xsh[npoints + i] = (*it)->GetX()[npoints-i-1];
	  ysh[npoints + i] = (*it)->GetY()[npoints-i-1] - (*it)->GetErrorYlow(npoints-i-1);
	}

      TGraph *centr = new TGraph(npoints, val_x, val_y);
      TGraph *high = new TGraph(npoints, val_x, val_high_y);
      TGraph *low = new TGraph(npoints, val_x, val_low_y);
      TGraph *shade = new TGraph(2*npoints, xsh, ysh);

      //Calculate maximum and minimum of y axis
      for (int i = 0; i < (*it)->GetN(); i++)
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

      //Set border lines and shade fill
      centr->SetLineColor((*it)->GetLineColor());
      centr->SetLineStyle(1);
      centr->SetLineWidth(opts.lwidth);
      high->SetLineColor((*it)->GetLineColor());
      high->SetLineStyle(1);
      high->SetLineWidth(opts.lwidth);
      low->SetLineColor((*it)->GetLineColor());
      low->SetLineStyle(1);
      low->SetLineWidth(opts.lwidth);
      shade->SetLineColor((*it)->GetLineColor());
      shade->SetFillColor((*it)->GetLineColor());
      shade->SetFillStyle((*it)->GetFillStyle());
      shade->SetLineWidth(0);

      //add graphs
      mg->Add((*it));
      mg_lines->Add(centr);
      mg_lines->Add(high);
      mg_lines->Add(low);
      mg_shade->Add(shade);
    }

  //Make the TCanvas
  TCanvas *cnv = new TCanvas(cnvname, "", opts.resolution, opts.resolution);
  cnvs.push_back(cnv);
  cnv->cd();
  if (opts.logx)
    cnv->SetLogx();
  cnv->SetLeftMargin(lmarg);
  cnv->SetRightMargin(rmarg);
  cnv->SetTopMargin(tmarg);

  //graphical settings
  mg->SetTitle(((string)" ; x  ; x" + pdflabels[ipdf] + "(x,Q^{2})").c_str());

  mg->Draw("AXIS"); //need to draw with A option to create axis

  //Set maximum and minimum
  if (mx != 0 || mn != 0)
    {
      double delta = mx - mn;
      mx = mx + delta * (0.1 + 0.1 * pdfgraphs.size());
      mn = mn - delta * 0.05;
      mg->SetMaximum(mx);
      mg->SetMinimum(mn);
    }
  
  mg->GetXaxis()->Set(100, opts.xmin, opts.xmax);
  mg->GetXaxis()->SetRange(opts.xmin, opts.xmax);
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

  //mg->Draw("LE3");
  mg_shade->Draw("f");
  mg_lines->Draw("L");

  //Make legend
  TLegend * leg = new TLegend(lmarg+0.05, 1-tmarg-0.05-pdfgraphs.size()*0.05, lmarg+0.5, 1-tmarg-0.01);
  leg->SetTextFont(62);
  leg->SetTextSize(txtsize);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  //  leg->AddEntry((TObject*)0, ((string)"x" + pdflabels[ipdf] + " - Q^{2} = " + q2str + " GeV^{2}").c_str(), "");
  leg->AddEntry((TObject*)0, ((string)"Q^{2} = " + q2str + " GeV^{2}").c_str(), "");

  for (vector <TGraphAsymmErrors*>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    {
      if (*it == 0)
	continue;
      if (opts.dobands)
	leg->AddEntry((*it), labels[it-pdfgraphs.begin()].c_str(), "lf");
      else
	leg->AddEntry((*it), labels[it-pdfgraphs.begin()].c_str(), "l");
    }

  leg->Draw();
  DrawLogo()->Draw();

  //--------------------------------------
  //Ratio Canvas
  //--------------------------------------
  if (opts.dirs.size() == 1)
    return cnvs;

  //Compute ratio graphs
  vector <TGraphAsymmErrors*> rlist;
  vector <TGraphAsymmErrors*>::iterator fit = pdfgraphs.begin();
  while (*fit == 0 && fit != pdfgraphs.end())
    fit++;
  for (vector <TGraphAsymmErrors*>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    {
      if (*it == 0)
	{
	  rlist.push_back(0);
	  continue;
	}

      TGraphAsymmErrors *r = (TGraphAsymmErrors*)(*it)->Clone();
      r->SetName(((string)(*it)->GetName() + "_ratio").c_str());
      for (int i = 0; i < (*it)->GetN(); i++)
	{
	  double ratio, rathigh, ratlow;

	  double val =  (*it)->GetY()[i];
	  double ref =  (*fit)->GetY()[i];
	  if (ref != 0)
	    ratio = val/ref;
	  else
	    ratio = 1;

	  if (opts.relerror)
	    ratio = 1;

	  if (opts.abserror)
	    ratio = 0;

	  double errhigh =  (*it)->GetErrorYhigh(i);
	  double errlow =  (*it)->GetErrorYlow(i);

	  rathigh = ( ref != 0)? (errhigh/ref) : 0;
	  ratlow = ( ref != 0)? (errlow/ref) : 0;
	  if (opts.relerror)
	    {
	      rathigh = ( val != 0 )? (errhigh/val) : 0;
	      ratlow = ( val != 0 )? (errlow/val) : 0;
	    }
	  if (opts.abserror)
	    {
	      rathigh = errhigh;
	      ratlow = errlow;
	    }
	  r->SetPoint(i, (*it)->GetX()[i], ratio);
	  r->SetPointError(i, 0, 0, ratlow, rathigh);
	}
      rlist.push_back(r);
    }

  mx = 1;
  mn = 1;
  if (opts.abserror)
    {
      mx = 0;
      mn = 0;
    }
  //Calculate maximum and minimum of y axis
  for (vector <TGraphAsymmErrors*>::iterator it = rlist.begin(); it != rlist.end(); it++)
    {
      if (*it == 0)
	continue;

      TGraphAsymmErrors* r = *it;
      double xmnforbds, xmxforbds;
      if (opts.logx)
	{
	  double xaxlength = opts.xmax / opts.xmin;
	  xmnforbds = opts.xmin * pow(xaxlength,1./4.);
	  xmxforbds = opts.xmax / pow(xaxlength,1./5.);
	}
      else
	{
	  double xaxlength = opts.xmax - opts.xmin;
	  xmnforbds = opts.xmin + xaxlength / 4.;
	  xmxforbds = opts.xmax - xaxlength / 5.;
	}
      if (opts.rmax == 0 && opts.rmin == 0)
	for (int i = 0; i < r->GetN(); i++)
	  {
	    double xi = r->GetX()[i];
	    double yi_h = r->GetY()[i] + r->GetErrorYhigh(i);
	    double yi_l = r->GetY()[i] - r->GetErrorYlow(i);
	    if (xi >= xmnforbds && xi <= xmxforbds)
	      {
		mx = max(mx, yi_h);
		mn = min(mn, yi_l);
	      }
	  }
      else
	{
	  mx = opts.rmax;
	  mn = opts.rmin;
	}
    }

  if (opts.rmax == 0 && opts.rmin == 0)
    if ((opts.abserror && (mx != 0 || mn != 0)) || (!opts.abserror && (mx != 1 || mn != 1)))
      {
	double delta = mx - mn;
	mx = mx + delta * (0.15 + 0.11 * pdfgraphs.size());
	mn = mn - delta * 0.3;
      }

  //prepare TGraphs for line borders and graph shade
  TMultiGraph * mg_ratio = new TMultiGraph(((string)cnvname + "_multigraph_ratio").c_str(), "");
  TMultiGraph * mg_ratio_lines = new TMultiGraph(((string)cnvname + "_multigraph_ratio_lines").c_str(), "");
  TMultiGraph * mg_ratio_shade = new TMultiGraph(((string)cnvname + "_multigraph_ratio_shade").c_str(), "");
  double tolerance = 0.01;   //tolerance for graph boundaries
  for (vector <TGraphAsymmErrors*>::iterator it = rlist.begin(); it != rlist.end(); it++)
    {
      if (*it == 0)
	continue;

      TGraphAsymmErrors* r = *it;

      int npoints = r->GetN();
      double val_x[npoints], val_y[npoints], val_high_y[npoints], val_low_y[npoints]; 
      double xsh[2*npoints], ysh[2*npoints];
      
      for (int i = 0; i < npoints; i++)
	{
	  //Set graphical safety boundaries
	  double val = r->GetY()[i];

	  double ratio = r->GetY()[i];
	  double high = r->GetY()[i] + r->GetErrorYhigh(i);
	  double low = r->GetY()[i] - r->GetErrorYlow(i);

	  double ratio_tol = ratio;
	  double high_tol = high;
	  double low_tol = low;

	  double delta = mx - mn;
	  if (ratio > (mx + delta * -tolerance))
	    ratio = mx + delta * -tolerance;
	  if (high > (mx + delta * -tolerance))
	    high = mx + delta * -tolerance;
	  if (low > (mx + delta * -tolerance))
	    low = mx + delta * -tolerance;

	  if (ratio_tol > (mx + delta * tolerance))
	    ratio_tol = mx + delta * tolerance;
	  if (high_tol > (mx + delta * tolerance))
	    high_tol = mx + delta * tolerance;
	  if (low_tol > (mx + delta * tolerance))
	    low_tol = mx + delta * tolerance;

	  if (ratio < (mn - delta * -tolerance))
	    ratio = mn - delta * -tolerance;
	  if (high < (mn - delta * -tolerance))
	    high = mn - delta * -tolerance ;
	  if (low < (mn - delta * -tolerance))
	    low = mn - delta * -tolerance;

	  if (ratio_tol < (mn - delta * tolerance))
	    ratio_tol = mn - delta * tolerance;
	  if (high_tol < (mn - delta * tolerance))
	    high_tol = mn - delta * tolerance;
	  if (low_tol < (mn - delta * tolerance))
	    low_tol = mn - delta * tolerance;

	  double errhigh = high - ratio;
	  double errlow = ratio - low;
	  r->SetPoint(i, r->GetX()[i], ratio);
	  r->SetPointError(i, 0, 0, errlow, errhigh);

	  val_x[i] = r->GetX()[i];
	  val_y[i] = ratio_tol;
	  val_high_y[i] = high_tol;
	  val_low_y[i] = low_tol;
	}
    
      //shade TGraph
      for (int i = 0; i < r->GetN(); i++)
	{
	  xsh[i] = r->GetX()[i];
	  ysh[i] = r->GetY()[i] + r->GetErrorYhigh(i);
	  xsh[npoints + i] = r->GetX()[npoints-i-1];
	  ysh[npoints + i] = r->GetY()[npoints-i-1] - r->GetErrorYlow(npoints-i-1);
	}

      TGraph *r_centr = new TGraph(npoints, val_x, val_y);
      TGraph *r_high = new TGraph(npoints, val_x, val_high_y);
      TGraph *r_low = new TGraph(npoints, val_x, val_low_y);
      TGraph *r_shade = new TGraph(2*npoints, xsh, ysh);

      //Set border lines and shade fill
      r_centr->SetLineColor(r->GetLineColor());
      r_centr->SetLineStyle(1);
      r_centr->SetLineWidth(opts.lwidth);
      r_high->SetLineColor(r->GetLineColor());
      r_high->SetLineStyle(1);
      r_high->SetLineWidth(opts.lwidth);
      r_low->SetLineColor(r->GetLineColor());
      r_low->SetLineStyle(1);
      r_low->SetLineWidth(opts.lwidth);
      r_shade->SetLineColor(r->GetLineColor());
      r_shade->SetLineColor(0);
      r_shade->SetFillColor(r->GetLineColor());
      r_shade->SetFillStyle(r->GetFillStyle());
      r_shade->SetLineWidth(0);

      //add graphs
      mg_ratio->Add(r);
      mg_ratio_lines->Add(r_centr);
      mg_ratio_lines->Add(r_high);
      mg_ratio_lines->Add(r_low);
      mg_ratio_shade->Add(r_shade);
    }

  //Make the TCanvas
  TCanvas *cnvr = new TCanvas(((string)cnvname + "_ratio").c_str(), "pdf", opts.resolution, opts.resolution);
  cnvs.push_back(cnvr);
  cnvr->cd();
  if (opts.logx)
    cnvr->SetLogx();
  cnvr->SetLeftMargin(lmarg);
  cnvr->SetRightMargin(rmarg);
  cnvr->SetTopMargin(tmarg);

  //graphical settings
  mg_ratio->Draw("AXIS"); //Create axis
  if ((opts.abserror && (mx != 0 || mn != 0)) || (!opts.abserror && (mx != 1 || mn != 1)))
    {
      mg_ratio->SetMaximum(mx);
      mg_ratio->SetMinimum(mn);
    }
  else
    {
      if (opts.abserror)
	{
	  mg_ratio->SetMaximum(1);
	  mg_ratio->SetMinimum(-1);
	}
      else
	{
	  mg_ratio->SetMaximum(2);
	  mg_ratio->SetMinimum(0);
	}
    }

  mg_ratio->GetXaxis()->SetTitle(" x  ");
  mg_ratio->GetYaxis()->SetTitle(((string)" x" + pdflabels[ipdf] + "(x,Q^{2})/x" + pdflabels[ipdf] + "(x,Q^{2})_{ref}").c_str());
  if (opts.relerror)
    mg_ratio->GetYaxis()->SetTitle(((string)" #deltax" + pdflabels[ipdf] + "/#deltax" + pdflabels[ipdf] + "_{ref}").c_str());
  if (opts.abserror)
    mg_ratio->GetYaxis()->SetTitle(((string)" #deltax" + pdflabels[ipdf] + "").c_str());


  mg_ratio->GetXaxis()->Set(100, opts.xmin, opts.xmax);
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
  mg_ratio->GetYaxis()->SetNdivisions(506);

  //  mg_ratio->Draw("ALE3");
  mg_ratio_shade->Draw("f");
  mg_ratio_lines->Draw("l");

  //Make legend
  TLegend * leg2 = new TLegend(lmarg+0.20, 1-tmarg-0.05-pdfgraphs.size()*0.05, lmarg+0.65, 1-tmarg-0.01);
  leg2->SetTextFont(62);
  leg2->SetTextSize(txtsize);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->AddEntry((TObject*)0, ((string)"Q^{2} = " + q2str + " GeV^{2}").c_str(), "");

  for (vector <TGraphAsymmErrors*>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    {
      if (*it == 0)
	continue;
      if (opts.dobands)
	leg2->AddEntry((*it), labels[it-pdfgraphs.begin()].c_str(), "lf");
      else
	leg2->AddEntry((*it), labels[it-pdfgraphs.begin()].c_str(), "l");
    }

  leg2->Draw();

  DrawLogo("dc")->Draw();

  return cnvs;
}

