#include <DataPainter.h>
#include <CommandParser.h>

#include <DrawLogo.h>

#include <TH1F.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TImage.h>

#include <iostream>
#include <math.h>

dataseth::dataseth(string dataname, string dir, string lab,
		   vector <float> bins1, vector <float> bins2, 
		   vector <float> data, vector <float> uncorerr, vector <float> toterr, 
		   vector <float> theory, vector <float> theoryshifted, 
		   vector <float> therrup, vector <float> therrdown, 
		   vector <float> pulls, bool Logx, bool Logy, float xmin, float xmax,
		   string xlabel, string ylabel)
  : name(dataname), label(lab), logx(Logx), logy(Logy)
{
  float bin[bins1.size() + 1];
  int i = 0;
  for (vector<float>::iterator it = bins1.begin(); it != bins1.end(); it++)
    {
      bin[i] = *it;
      i++;
    }
  bin[i] = *(bins2.end()-1);
  hdata = new TH1F((name + dir +"_data").c_str(), "", bins1.size(),  bin);
  hdatatot = new TH1F((name + dir + "_datatot").c_str(), "", bins1.size(),  bin);
  hth = new TH1F((name + dir + "_th").c_str(), "", bins1.size(),  bin);
  hthshift = new TH1F((name + dir + "_thshift").c_str(), "", bins1.size(),  bin);
  htherr = new TH1F((name + dir + "_therr").c_str(), "", bins1.size(),  bin);
  htherrup = new TH1F((name + dir + "_therrup").c_str(), "", bins1.size(),  bin);
  htherrdown = new TH1F((name + dir + "_therrdown").c_str(), "", bins1.size(),  bin);
  hpull = new TH1F((name + dir + "_pull").c_str(), "", bins1.size(),  bin);

  if (xmin != 0 && xmax != 0)
    {
      hdata->SetAxisRange(xmin, xmax);
      hdatatot->SetAxisRange(xmin, xmax);
      hth->SetAxisRange(xmin, xmax);
      hthshift->SetAxisRange(xmin, xmax);
      htherr->SetAxisRange(xmin, xmax);
      htherrup->SetAxisRange(xmin, xmax);
      htherrdown->SetAxisRange(xmin, xmax);
      hpull->SetAxisRange(xmin, xmax);
    }

  //  hdata->SetXTitle(xlabel.c_str());
  hdatatot->SetXTitle(xlabel.c_str());
  hpull->SetXTitle(xlabel.c_str());
  hdatatot->SetYTitle(ylabel.c_str());

  //set rapidity as default label 
  if (xlabel == "" && ylabel == "")
    {
      hdatatot->SetXTitle("y");
      hpull->SetXTitle("y");
      hdatatot->SetYTitle("d#sigma/dy [pb]");
    }

  for (unsigned int b = 0; b < data.size(); b++)
    {
      hdata->SetBinContent(b + 1, data[b]);
      hdata->SetBinError(b + 1, uncorerr[b]);
      hdatatot->SetBinContent(b + 1, data[b]);
      hdatatot->SetBinError(b + 1, toterr[b]);
      hth->SetBinContent(b + 1, theory[b]);
      hthshift->SetBinContent(b + 1, theoryshifted[b]);
      htherr->SetBinContent(b + 1, theory[b] + (therrup[b] - therrdown[b]) / 2);
      htherr->SetBinError(b + 1, (therrup[b] + therrdown[b]) / 2);
      htherrup->SetBinContent(b + 1, theory[b] + therrup[b]);
      htherrdown->SetBinContent(b + 1, theory[b] - therrdown[b]);
      //invert pulls -> (theory - data)
      hpull->SetBinContent(b + 1, -pulls[b]);
    }
}

TCanvas * DataPainter(int dataindex, vector <dataseth> datahistos)
{
  if (datahistos.size() < 1)
    {
      cout << "Empty dataset vector for dataset idx: " << dataindex << endl;
      return 0;
    }

  char cnvname[10];
  sprintf(cnvname, "%d_pulls",  dataindex);

  TCanvas * cnv = new TCanvas(cnvname, "", 0, 0, 2 * opts.resolution, opts.resolution);
  cnv->cd();

  TH1F * data = datahistos[0].getdata();
  TH1F * datatot = datahistos[0].getdatatot();
  string dataname = datahistos[0].getname();

  cnv->Divide(2, 1);
  cnv->GetPad(1)->SetLeftMargin(lmarg);
  if (datahistos[0].getlogy())
    cnv->GetPad(1)->SetLogy();
  if (datahistos[0].getlogx())
    cnv->GetPad(1)->SetLogx();
  cnv->cd(1);

  datatot->GetYaxis()->SetLabelFont(62);
  datatot->GetYaxis()->SetTitleFont(62);
  datatot->GetYaxis()->SetLabelSize(txtsize);
  datatot->GetYaxis()->SetTitleSize(txtsize);
  datatot->GetYaxis()->SetTitleOffset(offset);

  datatot->GetXaxis()->SetLabelFont(62);
  datatot->GetXaxis()->SetTitleFont(62);
  datatot->GetXaxis()->SetLabelSize(txtsize);
  datatot->GetXaxis()->SetTitleSize(txtsize);

  datatot->SetStats(0);
  datatot->SetFillColor(kYellow);

  //Evaluate maximum and minimum
  TH1F * dataerr = (TH1F*) datatot->Clone();
  for (int b = 1; b <= datatot->GetNbinsX(); b++)
    dataerr->SetBinContent(b, datatot->GetBinContent(b) + datatot->GetBinError(b));
  float mx = dataerr->GetMaximum();
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      mx = max(mx, (float)((*it).gettherrup()->GetMaximum()));
      mx = max(mx, (float)((*it).getthshift()->GetMaximum()));
    }
  for (int b = 1; b <= datatot->GetNbinsX(); b++)
    dataerr->SetBinContent(b, datatot->GetBinContent(b) - datatot->GetBinError(b));
  float mn = dataerr->GetMinimum();
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      mn = min(mn, (float)((*it).gettherrdown()->GetMinimum()));
      mn = min(mn, (float)((*it).getthshift()->GetMinimum()));
    }

  if (datahistos[0].getlogy())
    {
      if (mn < 0)
	mn = 0.000001;
      float ratio = mx / mn;
      mx = mx * pow(10, log10(ratio) * 0.35);
      mn = mn / pow(10, log10(ratio) * 0.65);
    }
  else
    {
      float delta = mx - mn;
      mx = mx + delta * 0.35;
      mn = mn - delta * 0.65;
    }

  datatot->SetMaximum(mx);
  datatot->SetMinimum(mn);

  datatot->Draw("e3");
  data->SetLineColor(1);
  data->SetMarkerStyle(20);
  data->SetMarkerSize(2 * opts.resolution / 1200);
  data->Draw("e1 same");

  //Main legend
  TLegend * leg = new TLegend(0.17, 0.13, 0.6, 0.35);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextAlign(12);
  leg->SetTextSize(txtsize * 0.8);
  leg->SetTextFont(62);
  leg->AddEntry(data, dataname.c_str(), "pl");
  leg->AddEntry(data, "#sigma uncorrelated", "pe");
  leg->AddEntry(datatot, "#sigma total", "f");
  TLine *cont = new TLine(0, 1, 1, 1);
  cont->SetLineStyle(1);
  TLine *dash = new TLine(0, 1, 1, 1);
  dash->SetLineStyle(2);
  leg->AddEntry(cont, "Theory", "l");
  leg->AddEntry(dash, "Theory + shifts", "l");

  //Auxiliary legend
  TLegend * leg2 = new TLegend(0.6, 0.13, 0.89, 0.13 + datahistos.size() * 0.045);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextAlign(12);
  leg2->SetTextFont(62);
  leg2->SetTextSize(txtsize * 0.8);

  //Plot theories
  int colindx = 0;
  int markindx = 0;
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      TGraphAsymmErrors * gtherr = new TGraphAsymmErrors((*it).getth());
      (*it).getthshift()->SetLineColor(opts.colors[colindx]);
      (*it).getthshift()->SetLineStyle(2);
      (*it).getthshift()->Draw("l same");

      (*it).getth()->SetLineColor(opts.colors[colindx]);
      if (!opts.points) //plot as continous line with dashed error bands
	{
	  (*it).getth()->Draw("l same");
	  if (opts.therr)
	    {    
	      (*it).gettherr()->SetLineColor(opts.colors[colindx]);
	      (*it).gettherr()->SetFillColor(opts.colors[colindx]);
	      (*it).gettherr()->SetFillStyle(opts.styles[colindx]);
	      float toterr = 0;
	      for (int b = 1; b <= (*it).gettherr()->GetNbinsX(); b++)
		toterr += (*it).gettherr()->GetBinError(b);
	      if (toterr > 0)
		(*it).gettherr()->Draw("e3 l same");
	    }
	}
      else //plot as tilted points with vertical error line
	{
	  gtherr->SetMarkerStyle(opts.markers[markindx]);
	  gtherr->SetLineColor(opts.colors[colindx]);
	  gtherr->SetMarkerSize(2 * opts.resolution / 1200);
	  gtherr->SetMarkerColor(opts.colors[colindx]);
	  for (int b = 0; b < gtherr->GetN(); b++)
	    {
	      //Set X error to 0
	      gtherr->SetPointEXlow(b, 0);
	      gtherr->SetPointEXhigh(b, 0);

	      //tilt horizontally
	      double x, y;
	      gtherr->GetPoint(b, x, y);
	      float width = (*it).getth()->GetBinWidth(b + 1);
	      float lowedge = (*it).getth()->GetBinLowEdge(b + 1);
	      x = lowedge + (it - datahistos.begin() + 1) * width/(datahistos.size() + 1);
	      gtherr->SetPoint(b, x, y);

	      //Set Y error
	      float errup, errdown;
	      if (opts.therr)
		{    
		  errup = (*it).gettherrup()->GetBinContent(b + 1) - (*it).getth()->GetBinContent(b + 1);
		  errdown = (*it).getth()->GetBinContent(b + 1) - (*it).gettherrdown()->GetBinContent(b + 1);
		}
	      else
		{
		  errup = 0;
		  errdown = 0;
		}
	      gtherr->SetPointEYhigh(b, errup);
	      gtherr->SetPointEYlow(b, errdown);
	    }
	  gtherr->Draw("P same");
	}
      colindx++;
      markindx++;
      if (!opts.points)
	if (opts.therr)
	  leg2->AddEntry((*it).gettherr(), ((*it).getlabel()).c_str(), "lf");
	else
	  leg2->AddEntry((*it).getth(), ((*it).getlabel()).c_str(), "l");
      else
	if (opts.therr)
	  leg2->AddEntry(gtherr, ((*it).getlabel()).c_str(), "pe");
	else
	  leg2->AddEntry(gtherr, ((*it).getlabel()).c_str(), "p");
    }

  //draw theory error borders
  colindx = 0;
  if (opts.therr)
    for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
      {
	(*it).gettherrup()->SetLineColor(opts.colors[colindx]);
	(*it).gettherrdown()->SetLineColor(opts.colors[colindx]);
	if (!opts.points)
	  {
	    (*it).gettherrup()->Draw("l same");
	    (*it).gettherrdown()->Draw("l same");
	  }
	colindx++;
      }
  leg->Draw();
  leg2->Draw();
  DrawLogo()->Draw();

  //Theory/Data ratio Pad
  cnv->GetPad(2)->Divide(1, 2);
  cnv->GetPad(2)->GetPad(1)->SetPad(0, 0.5, 1, 1);
  cnv->GetPad(2)->GetPad(1)->SetTopMargin(0.2);
  cnv->GetPad(2)->GetPad(1)->SetLeftMargin(lmarg);
  cnv->GetPad(2)->GetPad(1)->SetBottomMargin(0.);
  if (datahistos[0].getlogx())
    cnv->GetPad(2)->GetPad(1)->SetLogx();

  cnv->GetPad(2)->cd(1);
  TH1F * refdata = (TH1F*)datahistos[0].getdata()->Clone();
  for (int b = 1; b <= refdata->GetNbinsX(); b++)
    refdata->SetBinError(b, 0);
  TH1F * r_data = (TH1F*)datahistos[0].getdata()->Clone();
  TH1F * r_datatot = (TH1F*)datahistos[0].getdatatot()->Clone();
  r_data->Divide(refdata);
  r_datatot->Divide(refdata);

  r_datatot->GetYaxis()->SetLabelFont(62);
  r_datatot->GetYaxis()->SetTitleFont(62);
  r_datatot->GetYaxis()->SetLabelSize(txtsize * 2);
  r_datatot->GetYaxis()->SetTitleSize(txtsize * 2);
  r_datatot->GetYaxis()->SetTitleOffset(offset / 2);
  r_datatot->SetYTitle("Theory/Data   ");
  r_datatot->GetYaxis()->SetNdivisions(504);

  r_datatot->GetXaxis()->SetLabelFont(62);
  r_datatot->GetXaxis()->SetTitleFont(62);
  r_datatot->GetXaxis()->SetLabelSize(txtsize * 2);
  r_datatot->GetXaxis()->SetTitleSize(txtsize * 2);
  //  r_datatot->GetXaxis()->SetTitleOffset(1);

  r_datatot->SetStats(0);

  //Evaluate maximum and minimum
  TH1F * r_dataerr = (TH1F*) r_datatot->Clone();
  for (int b = 1; b <= r_datatot->GetNbinsX(); b++)
    r_dataerr->SetBinContent(b, r_datatot->GetBinContent(b) + r_datatot->GetBinError(b));
  mx = r_dataerr->GetBinContent(r_dataerr->GetMaximumBin());
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      TH1F * r_therrup = (TH1F*)(*it).gettherrup()->Clone();
      r_therrup->Divide(refdata);
      mx = max(mx, (float)(r_therrup->GetMaximum()));
      TH1F * r_thshift = (TH1F*)(*it).getthshift()->Clone();
      r_thshift->Divide(refdata);
      mx = max(mx, (float)(r_thshift->GetMaximum()));
    }
  for (int b = 1; b <= r_dataerr->GetNbinsX(); b++)
    r_dataerr->SetBinContent(b, r_datatot->GetBinContent(b) - r_datatot->GetBinError(b));
  mn = r_dataerr->GetBinContent(r_dataerr->GetMinimumBin());
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      TH1F * r_therrdown = (TH1F*)(*it).gettherrdown()->Clone();
      r_therrdown->Divide(refdata);
      mn = min(mn, (float)(r_therrdown->GetMinimum()));
      TH1F * r_thshift = (TH1F*)(*it).getthshift()->Clone();
      r_thshift->Divide(refdata);
      mn = min(mn, (float)(r_thshift->GetMinimum()));
    }
  float delta = mx - mn;
  r_datatot->SetMaximum(mx + delta * 0.2);
  r_datatot->SetMinimum(mn - delta * 0.2);

  //  r_datatot->SetMaximum(1.1);
  //  r_datatot->SetMinimum(0.9);

  //plot data
  r_datatot->Draw("e3");
  r_data->Draw("e1 same");

  //plot lines at 1
  TLine *r_one = new TLine(r_datatot->GetBinLowEdge(r_datatot->GetXaxis()->GetFirst()), 1, r_datatot->GetXaxis()->GetBinUpEdge(r_datatot->GetXaxis()->GetLast()), 1);
  r_one->SetLineStyle(2);
  r_one->SetLineStyle(1);
  r_one->Draw();

  //plot Ratios
  colindx = 0;
  markindx = 0;
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      //Prepare the ratio histograms
      TH1F * r_thshift = (TH1F*)(*it).getthshift()->Clone();
      for (int b = 1; b <= r_thshift->GetNbinsX(); b++)
      	r_thshift->SetBinError(b, 0);
      r_thshift->Divide(refdata);
      TH1F * r_th = (TH1F*)(*it).getth()->Clone();
      r_th->Divide(refdata);
      for (int b = 1; b <= r_th->GetNbinsX(); b++)
      	r_th->SetBinError(b, 0);
      TH1F * r_therr = (TH1F*)(*it).gettherr()->Clone();
      r_therr->Divide(refdata);
      TH1F * r_therrup = (TH1F*)(*it).gettherrup()->Clone();
      r_therrup->Divide(refdata);
      TH1F * r_therrdown = (TH1F*)(*it).gettherrdown()->Clone();
      r_therrdown->Divide(refdata);
      for (int b = 1; b <= r_therr->GetNbinsX(); b++)
	r_therr->SetBinError(b, (r_therrup->GetBinContent(b) - r_therrdown->GetBinContent(b)) / 2 );

      //Draw
      r_thshift->Draw("l same");
      if (!opts.points) //plot as continous line with dashed error bands
	{
	  r_th->Draw("l same");
	  if (opts.therr)
	    {
	      float toterr = 0;
	      for (int b = 1; b <= (*it).gettherr()->GetNbinsX(); b++)
		toterr += (*it).gettherr()->GetBinError(b);
	      if (toterr > 0)
		r_therr->Draw("e3 l same");
	    }
	}
      else //plot as tilted TGraphs
	{
	  TGraphAsymmErrors * r_gtherr = new TGraphAsymmErrors(r_th);
	  r_gtherr->SetMarkerStyle(opts.markers[markindx]);
	  r_gtherr->SetLineColor(opts.colors[colindx]);
	  r_gtherr->SetMarkerSize(2 * opts.resolution / 1200);
	  r_gtherr->SetMarkerColor(opts.colors[colindx]);
	  for (int b = 0; b < r_gtherr->GetN(); b++)
	    {
	      //Set X error to 0
	      r_gtherr->SetPointEXlow(b, 0);
	      r_gtherr->SetPointEXhigh(b, 0);

	      //tilt horizontally
	      double x, y;
	      r_gtherr->GetPoint(b, x, y);
	      float width = r_th->GetBinWidth(b + 1);
	      float lowedge = r_th->GetBinLowEdge(b + 1);
	      x = lowedge + (it - datahistos.begin() + 1) * width/(datahistos.size() + 1);
	      r_gtherr->SetPoint(b, x, y);
	      //Set Y error
	      float errup, errdown;
	      if (opts.therr)
		{    
		  errup = r_therrup->GetBinContent(b + 1) - r_th->GetBinContent(b + 1);
		  errdown = r_th->GetBinContent(b + 1) - r_therrdown->GetBinContent(b + 1);
		}
	      else
		{
		  errup = 0;
		  errdown = 0;
		}
	      r_gtherr->SetPointEYhigh(b, errup);
	      r_gtherr->SetPointEYlow(b, errdown);
	    }
	  r_gtherr->Draw("P same");
	}
      colindx++;
      markindx++;
    }	  
  
  //draw theory error borders
  if (opts.therr)
    for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
      {
	TH1F * r_therrup = (TH1F*)(*it).gettherrup()->Clone();
	r_therrup->Divide(refdata);
	for (int b = 1; b <= r_therrup->GetNbinsX(); b++)
	  r_therrup->SetBinError(b, 0);

	TH1F * r_therrdown = (TH1F*)(*it).gettherrdown()->Clone();
	r_therrdown->Divide(refdata);
	for (int b = 1; b <= r_therrdown->GetNbinsX(); b++)
	  r_therrdown->SetBinError(b, 0);

	if (!opts.points)
	  {
	    r_therrup->Draw("l same");
	    r_therrdown->Draw("l same");
	  }
      }

  //Theory-Data pulls pad
  cnv->GetPad(2)->GetPad(2)->SetPad(0, 0, 1, 0.5);
  cnv->GetPad(2)->GetPad(2)->SetTopMargin(0.0);
  cnv->GetPad(2)->GetPad(2)->SetLeftMargin(lmarg);
  cnv->GetPad(2)->GetPad(2)->SetBottomMargin(0.2);
  if (datahistos[0].getlogx())
    cnv->GetPad(2)->GetPad(2)->SetLogx();

  cnv->GetPad(2)->cd(2);

  TH1F * pull = datahistos[0].getpull();

  pull->GetYaxis()->SetLabelFont(62);
  pull->GetYaxis()->SetTitleFont(62);
  pull->GetYaxis()->SetLabelSize(txtsize * 2);
  pull->GetYaxis()->SetTitleSize(txtsize * 2);
  pull->GetYaxis()->SetTitleOffset(offset / 2);
  pull->GetYaxis()->SetNdivisions(504);
  pull->SetYTitle("#frac{Theory+shifts - Data}{#sigma uncor}");

  pull->GetXaxis()->SetLabelFont(62);
  pull->GetXaxis()->SetTitleFont(62);
  pull->GetXaxis()->SetLabelSize(txtsize * 2);
  pull->GetXaxis()->SetTitleSize(txtsize * 2);
  //  pull->GetXaxis()->SetTitleOffset(1);

  pull->SetStats(0);
  pull->SetMinimum(-3.5);
  pull->SetMaximum(3.5);

  //plot axis
  pull->Draw("][");

  //plot lines at 1, -1, 0
  TLine *one = new TLine(pull->GetBinLowEdge(pull->GetXaxis()->GetFirst()), 1, pull->GetXaxis()->GetBinUpEdge(pull->GetXaxis()->GetLast()), 1);
  one->SetLineStyle(2);
  TLine *minusone = new TLine(pull->GetBinLowEdge(pull->GetXaxis()->GetFirst()), -1, pull->GetXaxis()->GetBinUpEdge(pull->GetXaxis()->GetLast()), -1);
  minusone->SetLineStyle(2);
  TLine *zero = new TLine(pull->GetBinLowEdge(pull->GetXaxis()->GetFirst()), 0, pull->GetXaxis()->GetBinUpEdge(pull->GetXaxis()->GetLast()), 0);
  one->Draw();
  minusone->Draw();
  zero->Draw();

  //plot pulls
  colindx = 0;
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      (*it).getpull()->SetFillColor(opts.colors[colindx]);
      if (datahistos.size() == 1)
	(*it).getpull()->SetFillStyle(1001);
      else
	(*it).getpull()->SetFillStyle(opts.styles[colindx]);
      (*it).getpull()->SetLineStyle(1);
      (*it).getpull()->SetLineWidth(2);
      (*it).getpull()->SetLineColor(opts.colors[colindx]);
      colindx++;
      (*it).getpull()->Draw("same ][");
    }	  
  //redraw lines over fill area
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      TH1F * redrawpull = (TH1F*)(*it).getpull()->Clone();
      redrawpull->SetFillStyle(0);
      redrawpull->Draw("same ][");
    }

  return cnv;
}
