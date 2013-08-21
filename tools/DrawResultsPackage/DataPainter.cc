#include <DataPainter.h>

#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>

#include <iostream>

dataseth::dataseth(string dataname, string dir,
		   vector <float> bins1, vector <float> bins2, 
		   vector <float> data, vector <float> uncorerr, vector <float> toterr, 
		   vector <float> theory, vector <float> theoryshifted, 
		   vector <float> pulls) 
  : name(dataname)
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

  for (unsigned int b = 0; b < data.size(); b++)
    {
      hdata->SetBinContent(b + 1, data[b]);
      hdata->SetBinError(b + 1, uncorerr[b]);
      hdatatot->SetBinContent(b + 1, data[b]);
      hdatatot->SetBinError(b + 1, toterr[b]);
      hth->SetBinContent(b + 1, theory[b]);
      hthshift->SetBinContent(b + 1, theoryshifted[b]);
      hpull->SetBinContent(b + 1, pulls[b]);
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

  TCanvas * cnv = new TCanvas(cnvname, "", 0, 0, 600, 600);
  cnv->cd();
  TH1F * data = datahistos[0].getdata();
  TH1F * datatot = datahistos[0].getdatatot();
  TH1F * th = datahistos[0].getth();
  TH1F * thshift = datahistos[0].getthshift();
  TH1F * pull = datahistos[0].getpull();
  string dataname = datahistos[0].getname();

  cnv->Divide(1, 2);
  cnv->GetPad(1)->SetPad(0, 0.33, 1, 1);
  cnv->GetPad(1)->SetBottomMargin(0);
  //  cnv->GetPad(1)->SetLogy();
  //  cnv->GetPad(1)->SetLogx();
  cnv->cd(1);
  datatot->GetYaxis()->SetLabelFont(62);
  datatot->GetYaxis()->SetLabelSize(0.05);
  datatot->GetYaxis()->SetTitleFont(62);
  datatot->GetYaxis()->SetTitleSize(0.05);
  datatot->GetYaxis()->SetTitleOffset(1);
  // datatot->SetYTitle("d#sigma/dM_{tt} [pb]");
  datatot->SetYTitle("d#sigma/dy [pb]");
  datatot->SetStats(0);
  datatot->SetFillColor(kYellow);
  float mx = datatot->GetMaximum();
  mx = max(mx, (float)(th->GetMaximum()));
  float mn = datatot->GetMinimum();
  mn = min(mn, (float)(th->GetMinimum()));
  mx = mx + (mx - mn) * 0.3;
  mn = mn - (mx - mn) * 0.7;
  datatot->SetMaximum(mx);
  datatot->SetMinimum(mn);
  datatot->Draw("e3");
  data->SetLineColor(1);
  data->SetMarkerStyle(8);
  data->SetMarkerSize(0.5);
  data->Draw("e1 same");
  th->SetLineColor(kRed);
  th->Draw("l same");
  thshift->SetLineColor(kRed);
  thshift->SetLineStyle(2);
  thshift->Draw("l same");
  TLegend * leg = new TLegend(0.12, 0.03, 0.6, 0.25);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextAlign(12);
  leg->SetTextSize(0.05);
  leg->SetTextFont(62);
  leg->AddEntry(data, dataname.c_str(), "pl");
  leg->AddEntry(data, "uncorrelated uncertainty", "pe");
  leg->AddEntry(datatot, "total uncertainty", "f");
  leg->AddEntry(th, "fit result", "l");
  leg->AddEntry(thshift, "fit result + shift", "l");

  int colindx = 0;
  int colors[] = {kBlue, kGreen+4};
  for (vector <dataseth>::iterator it = (datahistos.begin() + 1); it != datahistos.end(); it++)
    {
      TH1F * nth = (*it).getth();
      TH1F * nthshift = (*it).getthshift();
      nth->SetLineColor(colors[colindx]);
      nth->Draw("l same");
      nthshift->SetLineColor(colors[colindx]);
      nthshift->SetLineStyle(2);
      nthshift->Draw("l same");
      colindx++;
      leg->AddEntry(nth, "fit result", "l");
      leg->AddEntry(nthshift, "fit result + shift", "l");
    }
  leg->Draw();


  cnv->GetPad(2)->SetPad(0, 0, 1, 0.33);
  cnv->GetPad(2)->SetTopMargin(0);
  cnv->GetPad(2)->SetBottomMargin(0.2);
  //  cnv->GetPad(2)->SetLogx();
  cnv->cd(2);
  pull->GetYaxis()->SetLabelFont(62);
  pull->GetYaxis()->SetLabelSize(0.1);
  pull->GetYaxis()->SetTitleFont(62);
  pull->GetYaxis()->SetTitleSize(0.1);
  pull->GetYaxis()->SetTitleOffset(0.5);
  pull->GetXaxis()->SetLabelFont(62);
  pull->GetXaxis()->SetLabelSize(0.1);
  pull->GetXaxis()->SetTitleFont(62);
  pull->GetXaxis()->SetTitleSize(0.1);
  pull->GetXaxis()->SetTitleOffset(1);
  pull->SetStats(0);
  pull->SetMinimum(-3.5);
  pull->SetMaximum(3.5);
  pull->SetFillColor(kGreen+4);
  pull->SetLineColor(kGreen+4);
  pull->SetYTitle("pulls");
  //  pull->SetXTitle("M_{tt}");
  pull->SetXTitle("y");
  TLine *one = new TLine(pull->GetBinLowEdge(1), 1, pull->GetXaxis()->GetBinUpEdge(pull->GetNbinsX()), 1);
  one->SetLineStyle(2);
  TLine *minusone = new TLine(pull->GetBinLowEdge(1), -1, pull->GetXaxis()->GetBinUpEdge(pull->GetNbinsX()), -1);
  minusone->SetLineStyle(2);
  pull->Draw();
  one->Draw();
  minusone->Draw();
  pull->Draw("same");

  return cnv;
}
