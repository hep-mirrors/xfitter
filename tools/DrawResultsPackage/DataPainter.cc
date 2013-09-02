#include <DataPainter.h>
#include <CommandParser.h>

#include <DrawLogo.h>

#include <TH1F.h>
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
  hpull = new TH1F((name + dir + "_pull").c_str(), "", bins1.size(),  bin);

  if (xmin != 0 && xmax != 0)
    {
      hdata->SetAxisRange(xmin, xmax);
      hdatatot->SetAxisRange(xmin, xmax);
      hth->SetAxisRange(xmin, xmax);
      hthshift->SetAxisRange(xmin, xmax);
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

  TCanvas * cnv = new TCanvas(cnvname, "", 0, 0, 4800, 2400);
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
      mx = max(mx, (float)((*it).getth()->GetMaximum()));
      mx = max(mx, (float)((*it).getthshift()->GetMaximum()));
    }
  for (int b = 1; b <= datatot->GetNbinsX(); b++)
    dataerr->SetBinContent(b, datatot->GetBinContent(b) - datatot->GetBinError(b));
  float mn = dataerr->GetMinimum();
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      mn = min(mn, (float)((*it).getth()->GetMinimum()));
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
  data->SetMarkerStyle(8);
  data->SetMarkerSize(3);
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
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      TH1F * nth = (*it).getth();
      TH1F * nthshift = (*it).getthshift();
      nth->SetLineColor(opts.colors[colindx]);
      nth->Draw("l same");
      nthshift->SetLineColor(opts.colors[colindx]);
      nthshift->SetLineStyle(2);
      nthshift->Draw("l same");
      colindx++;
      leg2->AddEntry(nth, ((*it).getlabel()).c_str(), "l");
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
      TH1F * r_th = (TH1F*)(*it).getth()->Clone();
      r_th->Divide(refdata);
      mx = max(mx, (float)(r_th->GetMaximum()));
      TH1F * r_thshift = (TH1F*)(*it).getthshift()->Clone();
      r_thshift->Divide(refdata);
      mx = max(mx, (float)(r_thshift->GetMaximum()));
    }
  for (int b = 1; b <= r_dataerr->GetNbinsX(); b++)
    r_dataerr->SetBinContent(b, r_datatot->GetBinContent(b) - r_datatot->GetBinError(b));
  mn = r_dataerr->GetBinContent(r_dataerr->GetMinimumBin());
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      TH1F * r_th = (TH1F*)(*it).getth()->Clone();
      r_th->Divide(refdata);
      mn = min(mn, (float)(r_th->GetMinimum()));
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
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      TH1F * r_th = (TH1F*)(*it).getth()->Clone();
      r_th->Divide(refdata);
      for (int b = 1; b <= r_th->GetNbinsX(); b++)
      	r_th->SetBinError(b, 0);
      r_th->Draw("l same");
      TH1F * r_thshift = (TH1F*)(*it).getthshift()->Clone();
      for (int b = 1; b <= r_thshift->GetNbinsX(); b++)
      	r_thshift->SetBinError(b, 0);
      r_thshift->Divide(refdata);
      r_thshift->Draw("l same");
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
