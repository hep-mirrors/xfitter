#include "DataPainter.h"
#include "CommandParser.h"
#include "DrawLogo.h"
#include "Outdir.h"

#include <TH1F.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>

#include <iostream>
#include <math.h>
#include <algorithm>

double hmin(TH1F *h)
{
  double min0 = h->GetBinContent(h->GetXaxis()->GetFirst()) + 1;
  for (int b = h->GetXaxis()->GetFirst(); b <= h->GetXaxis()->GetLast(); b++)
    if (h->GetBinContent(b) != 0)
      min0 = min(min0, h->GetBinContent(b));
  return min0;	
}

struct range
{
  double lowedge;
  double upedge;
};

vector <range> historanges(TH1F *h)
{
  vector <range> ranges;
  range temp;
  temp.lowedge = h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetFirst());
  int b = h->GetXaxis()->GetFirst();
  for (; b <= h->GetXaxis()->GetLast(); b++)
    if (h->GetBinContent(b) == 0)
      {
	temp.upedge = h->GetXaxis()->GetBinLowEdge(b - 1);
	ranges.push_back(temp);
	temp.lowedge = h->GetXaxis()->GetBinUpEdge(b);
      }
  temp.upedge = h->GetXaxis()->GetBinUpEdge(b - 2);
  ranges.push_back(temp);
  return ranges;
}

void Subplot::Draw(TH1F* histo, string opt)
{
  if (bins1.size() == 1)
    if (opt.find("E3") != string::npos)
      opt.replace(opt.find("E3"), opt.find("E3")+2, "E2");
      
  if (maketgraph)
    {
      TGraphAsymmErrors * graph = new TGraphAsymmErrors(histo);

      //Set correct x point
      vector <double> valy;
      for (int i = 0; i < graph->GetN(); i++)
	valy.push_back(graph->GetY()[i]);
      for (int i = 0; i < graph->GetN(); i++)
	{
	  graph->SetPoint(i, valx[i], valy[i]);
	  graph->SetPointEXhigh(i, 0);
	  graph->SetPointEXlow(i, 0);
	}
      graph->Sort();

      if (opt.find("same") == string::npos)
	opt.insert(0, "A");

      if (opt.find("E1") != string::npos)
	opt.erase(opt.find("E1")+1);

      if (opt.find("][") != string::npos) //this is a pull histo, force drawing as histogram
	{
	  vector <double> valy;
	  for (int i = 0; i < graph->GetN(); i++)
	    valy.push_back(graph->GetY()[i]);
	  for (int i = 0; i < graph->GetN(); i++)
	    histo->SetBinContent(i+1, valy[i]);
	  histo->Draw(opt.c_str());
	}
      else
	graph->Draw(opt.c_str());
    }
  else histo->Draw(opt.c_str());
}

TCanvas * DataPainter(int dataindex, int subplotindex)
{
  vector <Subplot> datahistos;
  vector <string> labels;
  for (vector<string>::iterator itl = opts.labels.begin(); itl != opts.labels.end(); itl++)
    if (datamap[*itl].datamap.find(dataindex) !=  datamap[*itl].datamap.end())
      if (datamap[*itl].datamap[dataindex].subplots.find(subplotindex) !=  datamap[*itl].datamap[dataindex].subplots.end())
	if (datamap[*itl].datamap[dataindex].subplots[subplotindex].IsValid())
	  {
	    datahistos.push_back(datamap[*itl].datamap[dataindex].subplots[subplotindex]);
	    labels.push_back(*itl);
	  }

  if (datahistos.size() < 1)
    return 0; //Empty dataset vector

  char cnvname[15];
  sprintf(cnvname, "data_%d-%d",  dataindex, subplotindex);

  if (opts.multitheory)
    {
      if (datahistos.size() == 2)
	{
	  opts.twopanels = true;
	  opts.threepanels = false;
	}
      if (datahistos.size() == 3)
	{
	  opts.twopanels = false;
	  opts.threepanels = true;
	}
      if (datahistos.size() > 3)
	{
	  cout << "Cannot plot in multitheory mode more than 3 directories" << endl;
	  return 0;
	}
    }
  
  TCanvas * cnv;
  if (opts.twopanels || opts.threepanels)
    cnv = new TCanvas(cnvname, "", 0, 0, 2 * opts.resolution, opts.resolution);
  else
    cnv = new TCanvas(cnvname, "", 0, 0, opts.resolution, opts.resolution);
  cnv->cd();

  TH1F * data = datahistos[0].getdata();
  TH1F * datatot = datahistos[0].getdatatot();
  char dtname[200];
  sprintf (dtname, "data_%d-%d", dataindex, subplotindex);
  string dataname = dtname;

  //Set the pads geometry
  float dy; //subpanels height

  if (opts.twopanels)
    dy = 0.5 * (1.-bmarg-tmarg);
  else if (opts.threepanels)
    dy = (1.-bmarg-tmarg)/3.;
  else //1 panel
    dy = (1.-bmarg-tmarg)/4.;
    

  TPad* Main;
  TPad* Ratio;
  TPad* Shifts;
  TPad* Pulls;
  float my, ry, sy, py; //my is the main panel height
  float mb;
  if (opts.twopanels || opts.threepanels)
    cnv->Divide(2, 1);
  else
    cnv->Divide(1, 2);

  Main = (TPad*)cnv->GetPad(1);

  if (opts.twopanels)
    cnv->GetPad(2)->Divide(1, 2);
  else if (opts.threepanels)
    cnv->GetPad(2)->Divide(1, 3);

  if (opts.twopanels || opts.threepanels)
    Ratio = (TPad*)cnv->GetPad(2)->GetPad(1);
  else
    Ratio = (TPad*)cnv->GetPad(2);


  //Main pad geometry
  Main->SetLeftMargin(lmarg+0.02);
  Main->SetRightMargin(rmarg);
  Main->SetTopMargin(tmarg);
  if (opts.twopanels || opts.threepanels)
    {
      Main->SetBottomMargin(bmarg);
      my = 1;
      mb = bmarg;
    }
  else
    {
      Main->SetPad(0, bmarg+dy, 1, 1);
      Main->SetTopMargin(tmarg/(1-dy-bmarg));
      Main->SetBottomMargin(marg0/(1-dy-bmarg));
      my = 1-dy-bmarg;
      mb = marg0/my;
    }

  
  //Ratio pad geometry
  if (opts.twopanels || opts.threepanels)
    {
      Ratio->SetPad(0, 1-tmarg-dy, 1, 1);
      Ratio->SetTopMargin(tmarg/(dy+tmarg));
      Ratio->SetBottomMargin(marg0/(dy+tmarg));
      ry = dy+tmarg;
    }
  else
    {
      Ratio->SetPad(0, 0, 1, bmarg+dy);
      Ratio->SetTopMargin(marg0/(bmarg+dy));
      Ratio->SetBottomMargin(bmarg/(bmarg+dy));
      ry = dy+bmarg;
    }
  Ratio->SetLeftMargin(lmarg+0.02);
  Ratio->SetRightMargin(rmarg);

  //Shifts pad geometry
  if (opts.threepanels)
    {
      Shifts = (TPad*)cnv->GetPad(2)->GetPad(2);
      Shifts->SetPad(0, bmarg+dy, 1, bmarg+dy+dy);
      Shifts->SetTopMargin(marg0/dy);
      Shifts->SetLeftMargin(lmarg+0.02);
      Shifts->SetRightMargin(rmarg);
      Shifts->SetBottomMargin(marg0/dy);
      sy = dy;
    }

  //Pulls pad geometry
  if (opts.twopanels || opts.threepanels)
    {
      int pullpad;
      if (opts.twopanels)
	pullpad = 2;
      if (opts.threepanels)
	pullpad = 3;

      Pulls = (TPad*)cnv->GetPad(2)->GetPad(pullpad);
      Pulls->SetPad(0, 0, 1, bmarg+dy);
      Pulls->SetTopMargin(marg0/(bmarg+dy));
      Pulls->SetLeftMargin(lmarg+0.02);
      Pulls->SetRightMargin(rmarg);
      Pulls->SetBottomMargin(bmarg/(bmarg+dy));
      py=dy+bmarg;
    }

  //Draw Main Pad
  Main->cd();
  if (datahistos[0].getlogy())
    Main->SetLogy();
  if (datahistos[0].getlogx())
    Main->SetLogx();

  //Set axis range to allow a small distance between border and plots
  float axmin = datahistos[0].getxmin();
  float axmax = datahistos[0].getxmax();
  int nbins = data->GetNbinsX();

  if (datahistos[0].getlogx())
    {
      float lenght = axmax / axmin;
      axmin = axmin / (pow(lenght, 1./20));
      axmax = axmax * (pow(lenght, 1./20));
    }
  else
    {
      float lenght = fabs(axmax - axmin);
      axmin = axmin - (lenght / 20);
      axmax = axmax + (lenght / 20);
    }

  //create template histogram for axis
  TH1F *up_templ = new TH1F(((string) "up_templ_" + cnvname).c_str(), "", nbins, axmin, axmax);

  up_templ->GetYaxis()->SetLabelFont(62);
  up_templ->GetYaxis()->SetTitleFont(62);
  up_templ->GetYaxis()->SetLabelSize(txtsize/my);
  up_templ->GetYaxis()->SetTitleSize(txtsize/my);
  up_templ->GetYaxis()->SetTitleOffset((offset+0.3) * my);
  up_templ->GetYaxis()->SetTitle(data->GetYaxis()->GetTitle());

  up_templ->GetXaxis()->SetLabelFont(62);
  up_templ->GetXaxis()->SetTitleFont(62);
  up_templ->GetXaxis()->SetLabelSize(txtsize/my);
  up_templ->GetXaxis()->SetTitleSize(txtsize/my);
  up_templ->GetXaxis()->SetTitle(data->GetXaxis()->GetTitle());

  up_templ->GetXaxis()->SetNdivisions(505);
  if (datahistos[0].getlogx())
    {
      up_templ->GetXaxis()->SetNoExponent();
      if (axmax/axmin < 90)
	up_templ->GetXaxis()->SetMoreLogLabels();
    }

  //Evaluate maximum and minimum
  float mx = 0;
  TH1F * dataerr = (TH1F*) datatot->Clone();
  for (int b = 1; b <= datatot->GetNbinsX(); b++)
    dataerr->SetBinContent(b, datatot->GetBinContent(b) + datatot->GetBinError(b));
  if (!opts.onlytheory)
    mx = dataerr->GetMaximum();
  for (vector <Subplot>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      if (opts.therr && !opts.noupband)
	mx = max(mx, (float)((*it).gettherrup()->GetMaximum()));
      else
	mx = max(mx, (float)((*it).getth()->GetMaximum()));
      if (!opts.onlytheory)
	mx = max(mx, (float)((*it).getthshift()->GetMaximum()));
    }

  float mn = mx;
  for (int b = 1; b <= datatot->GetNbinsX(); b++)
    dataerr->SetBinContent(b, datatot->GetBinContent(b) - datatot->GetBinError(b));
  if (!opts.onlytheory)
    mn = hmin(dataerr);
  for (vector <Subplot>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      if (opts.therr && !opts.noupband)
	mn = min(mn, (float)(hmin((*it).gettherrdown())));
      else
	mn= min(mn, (float)(hmin((*it).getth())));
      if (!opts.onlytheory)
	mn = min(mn, (float)(hmin((*it).getthshift())));
    }

  if (datahistos[0].getlogy())
    {
      if (mn < 0)
	mn = 0.000001;
      float ratio = mx / mn;
      mx = mx * pow(10, log10(ratio) * 0.45/my);
      mn = mn / pow(10, log10(ratio) * 0.7/my);
    }
  else
    {
      float delta = mx - mn;
      mx = mx + delta * 0.45/my;
      mn = mn - delta * 0.7/my;
    }
  if (datahistos[0].getymax() != 0 && datahistos[0].getymin() != 0)
    {
      mx = datahistos[0].getymax();
      mn = datahistos[0].getymin();
    }

  up_templ->SetMaximum(mx);
  up_templ->SetMinimum(mn);
  up_templ->SetStats(0);
  up_templ->Draw("AXIS");

  data->SetLineColor(1);
  data->SetMarkerStyle(20);
  data->SetMarkerSize(2 * opts.resolution / 1200);

  datatot->SetFillColor(opts.errbandcol);
  datatot->SetLineColor(opts.errbandcol);
  vector <range> dtranges = historanges(datatot);
  for (vector<range>::iterator r = dtranges.begin(); r != dtranges.end(); r++)
    {
      datatot->SetAxisRange((*r).lowedge, (*r).upedge);
      if (!opts.onlytheory)
	datahistos[0].Draw((TH1F*)datatot->Clone(), "PE3 same");
    }
  if (!opts.onlytheory)
    datahistos[0].Draw(data, "PE1 same");
  //reset axis range
  datatot->GetXaxis()->SetRange(datahistos[0].getlowrange(), datahistos[0].getuprange());

  if (datahistos[0].getextralabel() != "")
    {
      TLatex l;
      l.SetNDC();
      l.SetTextFont(42);

      float vertdist;
      float txtsz;
      string infolabel;
      if (opts.atlasinternal || opts.atlaspreliminary || opts.atlas)
	{
	  vertdist = 0.10;
	  txtsz = 1.;
	  infolabel = datahistos[0].getextralabel() + "; " + datahistos[0].getlumilabel();
	}
      else
	{
	  vertdist = 0.05;
	  txtsz = 1.;
	  infolabel = datahistos[0].getextralabel();
	}
      l.SetTextSize(txtsz*0.04/my);
      l.DrawLatex(lmarg+0.05, (1-tmarg/my) - vertdist/my, infolabel.c_str());
    }

  if (datahistos[0].getlumilabel() != "")
    if (!(opts.atlasinternal || opts.atlaspreliminary || opts.atlas))
      {
	TLatex l;
	l.SetNDC();
	l.SetTextFont(42);
	
	float vertdist;
	float txtsz;
	vertdist = 0.13;
	txtsz = 1.;
	l.SetTextSize(txtsz*0.04/my);
	l.DrawLatex(lmarg+0.05, (1-tmarg/my) - vertdist/my, datahistos[0].getlumilabel().c_str());
      }
  
  //Main legend
  TPaveText* leg1;
  if (opts.onlytheory)
    {
      leg1 = new TPaveText(lmarg+0.05, mb+0.03, lmarg+0.05+0.3, mb+0.03+0.04/my, "NDC");
      leg1->AddText(opts.theorylabel.c_str());
    }
  else
    {
      TLegend * leg;
      if (opts.nothshifts)
	leg = new TLegend(lmarg+0.04, mb+0.03, lmarg+0.04+0.30, mb+0.03+0.12/my);
      else
	leg = new TLegend(lmarg+0.04, mb+0.03, lmarg+0.04+0.30, mb+0.03+0.2/my);
      string datalab = (string) "Data " + datahistos[0].gettitle();
      if (datahistos[0].getexperiment() != "")
	datalab = datahistos[0].getexperiment() + " " + datalab;
      if (!opts.onlytheory)
	{
	  leg->AddEntry(data, datalab.c_str(), "pl");
	  leg->AddEntry(data, "#delta uncorrelated", "pe");
	  leg->AddEntry(datatot, "#delta total", "f");
	}
      TH1 *mark = (TH1F*)datahistos[0].getth()->Clone();
      mark->SetMarkerStyle(opts.markers[labels[0]]);
      mark->SetMarkerSize(2 * opts.resolution / 1200);
      mark->SetMarkerColor(kBlack);
      TLine *cont = new TLine(0, 1, 1, 1);
      cont->SetLineStyle(1);
      cont->SetLineWidth(opts.lwidth);
      TLine *dash = new TLine(0, 1, 1, 1);
      dash->SetLineStyle(2);
      dash->SetLineWidth(opts.lwidth);
      if (datahistos.size() == 1)
	{
	  cont->SetLineColor(opts.colors[labels[0]]);
	  dash->SetLineColor(opts.colors[labels[0]]);
	  mark->SetMarkerColor(opts.colors[labels[0]]);
	}
      if (opts.onlytheory)
	leg->AddEntry((TObject*)0, opts.theorylabel.c_str(), "");
      else
	{
	  if (opts.nothshifts)
	    if ((opts.points && !datahistos[0].bincenter()) || datahistos[0].nbins() == 1)
	      leg->AddEntry(mark, opts.theorylabel.c_str(), "p");
	    else
	      leg->AddEntry(cont, opts.theorylabel.c_str(), "l");
	  if (!opts.nothshifts)
	    leg->AddEntry(dash, (opts.theorylabel + " + shifts").c_str(), "l");
	}
      leg1 = (TPaveText*)leg;
    }
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextAlign(12);
  leg1->SetTextSize(txtsize * 0.8/my);
  leg1->SetTextFont(62);

  //Auxiliary legend
  int leg2size =  datahistos.size();
  if (leg2size == 1 && (opts.therr && !opts.noupband && (datahistos.begin())->HasTherr()))
    leg2size++;
  TLegend * leg2 = new TLegend(lmarg+0.4, mb+0.03, lmarg+0.7, mb+0.03 + leg2size * 0.04/my);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextAlign(12);
  leg2->SetTextFont(62);
  leg2->SetTextSize(txtsize * 0.8/my);

  //Plot theories
  for (vector <Subplot>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      TGraphAsymmErrors * gtherr = new TGraphAsymmErrors((*it).getth());
      (*it).getthshift()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
      (*it).getthshift()->SetLineStyle(2);
      (*it).getthshift()->SetLineWidth(opts.lwidth);

      vector <range> thranges = historanges((*it).getthshift());
      for (vector<range>::iterator r = thranges.begin(); r != thranges.end(); r++)
	{
	  (*it).getthshift()->SetAxisRange((*r).lowedge, (*r).upedge);
	  if (!opts.onlytheory)
	    if (!opts.nothshifts)
	      if ((!opts.points || (*it).bincenter()) && (*it).nbins() > 1) //plot as continous line
		(*it).Draw((TH1F*)(*it).getthshift()->Clone(), "LX same");
	      else
		(*it).Draw((TH1F*)(*it).getthshift()->Clone(), "hist ][ same"); //plot as histogram in points mode
	}
      (*it).getthshift()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());

      (*it).getth()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
      (*it).getth()->SetLineWidth(opts.lwidth);
      if (opts.bw)
	(*it).getth()->SetLineStyle(opts.lstyles[labels[it-datahistos.begin()]]);

      if ((!opts.points || (*it).bincenter()) && (*it).nbins() > 1) //plot as continous line with dashed error bands
	{
	  for (vector<range>::iterator r = thranges.begin(); r != thranges.end(); r++)
	    {
	      (*it).getth()->SetAxisRange((*r).lowedge, (*r).upedge);
	      (*it).Draw((TH1F*)(*it).getth()->Clone(), "LX same");
	    }
	  (*it).getth()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());

	  if (opts.therr && !opts.noupband)
	    {    
	      (*it).gettherr()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	      if (opts.bw)
		(*it).gettherr()->SetLineStyle(opts.lstyles[labels[it-datahistos.begin()]]);
	      (*it).gettherr()->SetMarkerSize(0);
	      (*it).gettherr()->SetFillColor(opts.colors[labels[it-datahistos.begin()]]);
	      (*it).gettherr()->SetFillStyle(opts.styles[labels[it-datahistos.begin()]]);
	      float toterr = 0;
	      for (int b = 1; b <= (*it).gettherr()->GetNbinsX(); b++)
		toterr += (*it).gettherr()->GetBinError(b);
	      if (toterr > 0)
		{
		  for (vector<range>::iterator r = thranges.begin(); r != thranges.end(); r++)
		    {
		      (*it).gettherr()->SetAxisRange((*r).lowedge, (*r).upedge);
		      (*it).Draw((TH1F*)(*it).gettherr()->Clone(), "E3L same");
		    }
		  (*it).gettherr()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
		}
	    }
	}
      else //plot as displaced points with vertical error line
	{
	  gtherr->SetMarkerStyle(opts.markers[labels[it-datahistos.begin()]]);
	  gtherr->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	  gtherr->SetMarkerSize(2 * opts.resolution / 1200);
	  gtherr->SetMarkerColor(opts.colors[labels[it-datahistos.begin()]]);
	  for (int b = 0; b < gtherr->GetN(); b++)
	    {
	      //Set X error to 0
	      gtherr->SetPointEXlow(b, 0);
	      gtherr->SetPointEXhigh(b, 0);

	      //displace horizontally
	      double x, y;
	      gtherr->GetPoint(b, x, y);
	      float width = (*it).getth()->GetBinWidth(b + 1);
	      float lowedge = (*it).getth()->GetBinLowEdge(b + 1);
	      x = lowedge + (it - datahistos.begin() + 1) * width/(datahistos.size() + 1);
	      gtherr->SetPoint(b, x, y);

	      //Set Y error
	      float errup, errdown;
	      if (opts.therr && !opts.noupband)
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
      if (datahistos.size() == 1)
	{
	  leg2->AddEntry((TObject*)0, (labels[it-datahistos.begin()]).c_str(), "");
	  if (opts.therr && !opts.noupband && (*it).HasTherr())
	    if ((!opts.points || (*it).bincenter()) && (*it).nbins() > 1)
	      leg2->AddEntry((*it).gettherr(), "Theory uncertainty", "lf");
	    else
	      leg2->AddEntry((*it).gettherr(), "Theory uncertainty", "pe");
	}
      else
	if ((!opts.points || (*it).bincenter()) && (*it).nbins() > 1)
	  if (opts.therr && !opts.noupband && (*it).HasTherr())
	    leg2->AddEntry((*it).gettherr(), (labels[it-datahistos.begin()]).c_str(), "lf");
	  else
	    leg2->AddEntry((*it).getth(), (labels[it-datahistos.begin()]).c_str(), "l");
	else
	  if (opts.therr && !opts.noupband && (*it).HasTherr())
	    leg2->AddEntry(gtherr, (labels[it-datahistos.begin()]).c_str(), "pe");
	  else
	    leg2->AddEntry(gtherr, (labels[it-datahistos.begin()]).c_str(), "p");
    }

  //draw theory error borders
  if (opts.therr && !opts.noupband)
    for (vector <Subplot>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
      {
	(*it).gettherrup()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	(*it).gettherrdown()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	(*it).gettherrup()->SetLineWidth(opts.lwidth);
	(*it).gettherrdown()->SetLineWidth(opts.lwidth);
	if ((!opts.points || (*it).bincenter()) && (*it).nbins() > 1)
	  {
	    vector <range> thranges = historanges((*it).getth());
	    for (vector<range>::iterator r = thranges.begin(); r != thranges.end(); r++)
	      {
		(*it).gettherrup()->SetAxisRange((*r).lowedge, (*r).upedge);
		(*it).Draw((TH1F*)(*it).gettherrup()->Clone(), "LX same");
		(*it).gettherrdown()->SetAxisRange((*r).lowedge, (*r).upedge);
		(*it).Draw((TH1F*)(*it).gettherrdown()->Clone(), "LX same");
	      }
	    (*it).gettherrup()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
	    (*it).gettherrdown()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
	  }
      }

  data->SetStats(0);
  data->SetLineColor(1);
  data->SetMarkerStyle(20);
  data->SetMarkerSize(2 * opts.resolution / 1200);
  if (!opts.onlytheory)
    datahistos[0].Draw(data, "PE1 same");

  leg1->Draw();
  leg2->Draw();

  //Theory/Data ratio Pad
  Ratio->cd();
  if (datahistos[0].getlogx())
    Ratio->SetLogx();

  TH1F * refdata = (TH1F*)datahistos[0].getref()->Clone();
  for (int b = 1; b <= refdata->GetNbinsX(); b++)
    refdata->SetBinError(b, 0);
  TH1F * r_data = (TH1F*)datahistos[0].getdata()->Clone();
  TH1F * r_datatot = (TH1F*)datahistos[0].getdatatot()->Clone();

  //Double counting of errors avoided using refdata which has 0 errors
  if (opts.diff)
    {
      r_data->Add(refdata, -1);
      r_datatot->Add(refdata, -1);
    }
  else
    {
      r_data->Divide(refdata);
      r_datatot->Divide(refdata);
    }

  //Ratio of theory over reference
  for (vector <Subplot>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      if (opts.diff)
	{
	  (*it).getrth()->Add(refdata, -1);
	  (*it).getrthshift()->Add(refdata, -1);
	  (*it).getrtherr()->Add(refdata, -1);
	  (*it).getrtherrup()->Add(refdata, -1);
	  (*it).getrtherrdown()->Add(refdata, -1);
	}
      else
	{
	  if (opts.onlytheory && opts.threlerr)
	    {
	      (*it).getrth()->Divide((*it).getth());
	      (*it).getrthshift()->Divide((*it).getth());
	      (*it).getrtherr()->Divide((*it).getth());
	      (*it).getrtherrup()->Divide((*it).getth());
	      (*it).getrtherrdown()->Divide((*it).getth());
	    }
	  else
	    {
	      (*it).getrth()->Divide(refdata);
	      (*it).getrthshift()->Divide(refdata);
	      (*it).getrtherr()->Divide(refdata);
	      (*it).getrtherrup()->Divide(refdata);
	      (*it).getrtherrdown()->Divide(refdata);
	    }
	}

      for (int b = 1; b <= (*it).getrth()->GetNbinsX(); b++)
	(*it).getrth()->SetBinError(b, 0);
      for (int b = 1; b <= (*it).getrthshift()->GetNbinsX(); b++)
	(*it).getrthshift()->SetBinError(b, 0);
      for (int b = 1; b <= (*it).getrtherr()->GetNbinsX(); b++)
	(*it).getrtherr()->SetBinError(b, ((*it).getrtherrup()->GetBinContent(b) - (*it).getrtherrdown()->GetBinContent(b)) / 2 );
      for (int b = 1; b <= (*it).getrtherrup()->GetNbinsX(); b++)
	(*it).getrtherrup()->SetBinError(b, 0);
      for (int b = 1; b <= (*it).getrtherrdown()->GetNbinsX(); b++)
	(*it).getrtherrdown()->SetBinError(b, 0);
    }

  //create template histogram for axis
  TH1F *r_templ = new TH1F(((string) "r_templ_" + cnvname).c_str(), "", nbins, axmin, axmax);

  r_templ->GetYaxis()->SetLabelFont(62);
  r_templ->GetYaxis()->SetTitleFont(62);
  r_templ->GetXaxis()->SetLabelFont(62);
  r_templ->GetXaxis()->SetTitleFont(62);
  r_templ->GetYaxis()->SetNdivisions(505);
  r_templ->GetXaxis()->SetNdivisions(505);
  if (datahistos[0].getlogx())
    {
      r_templ->GetXaxis()->SetNoExponent();
      if (axmax/axmin < 90)
	r_templ->GetXaxis()->SetMoreLogLabels();
    }

  r_templ->GetYaxis()->SetLabelSize(txtsize/ry);
  r_templ->GetYaxis()->SetTitleSize(txtsize/ry);
  r_templ->GetYaxis()->SetTitleOffset((offset+0.3) * ry);
  string ytitle = "";
  if (opts.diff)
    {
      if (opts.onlytheory)
	ytitle = "Difference";
      else if (opts.ratiototheory)
	ytitle = (string) "Data-" + opts.theorylabel;
      else
	ytitle = "Theory-Data";
    }
  else
    {
      if (opts.onlytheory)
	ytitle = "Ratio";
      else if (opts.ratiototheory)
	ytitle = (string) "Data/" + opts.theorylabel;
      else
	ytitle = "Theory/Data";
    }
  r_templ->SetYTitle(ytitle.c_str());

  r_templ->GetXaxis()->SetLabelSize(txtsize/ry);
  r_templ->GetXaxis()->SetTitleSize(txtsize/ry);
  //  r_templ->GetXaxis()->SetTitleOffset((offset+0.3) * ry);
  if (opts.twopanels || opts.threepanels)
    {
      r_templ->GetXaxis()->SetLabelSize(0);
      r_templ->GetXaxis()->SetTitleSize(0);
    }
  r_templ->GetXaxis()->SetTitle(data->GetXaxis()->GetTitle());

  r_templ->SetStats(0);

  //Evaluate maximum and minimum
  mx = 0;
  TH1F * r_dataerr = (TH1F*) r_datatot->Clone();
  for (int b = 1; b <= r_datatot->GetNbinsX(); b++)
    r_dataerr->SetBinContent(b, r_datatot->GetBinContent(b) + r_datatot->GetBinError(b));
  if (!opts.onlytheory)
    mx = r_dataerr->GetBinContent(r_dataerr->GetMaximumBin());
  for (vector <Subplot>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      if (opts.therr)
	{
	  mx = max(mx, (float)((*it).getrtherrup()->GetMaximum()));
	  mx = max(mx, (float)((*it).getrtherrdown()->GetMaximum()));
	}
      else
	mx = max(mx, (float)((*it).getrth()->GetMaximum()));
      if (!opts.threepanels)
	if (!opts.onlytheory)
	  mx = max(mx, (float)((*it).getrthshift()->GetMaximum()));
    }

  mn = mx;
  for (int b = 1; b <= r_dataerr->GetNbinsX(); b++)
    r_dataerr->SetBinContent(b, r_datatot->GetBinContent(b) - r_datatot->GetBinError(b));
  if (!opts.onlytheory)
    mn = hmin(r_dataerr);
  for (vector <Subplot>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      if (opts.therr)
	{
	  mn = min(mn, (float)(hmin((*it).getrtherrdown())));
	  mn = min(mn, (float)(hmin((*it).getrtherrup())));
	}
      else
	mn = min(mn, (float)(hmin((*it).getrth())));
      if (!opts.threepanels)
	if (!opts.onlytheory)
	  mn = min(mn, (float)(hmin((*it).getrthshift())));
    }
  float delta = mx - mn;
  if (datahistos[0].getymaxr() != 0)
    {
      mx = datahistos[0].getymaxr();
      mn = datahistos[0].getyminr();
      delta = 0;
    }

  r_templ->SetMaximum(mx + delta * 0.2);
  r_templ->SetMinimum(mn - delta * 0.2);

  //draw axis
  r_templ->DrawCopy("AXIS");

  //draw data
  vector <range> rdtranges = historanges(r_datatot);
  for (vector<range>::iterator r = rdtranges.begin(); r != rdtranges.end(); r++)
    {
      r_datatot->SetAxisRange((*r).lowedge, (*r).upedge);
      if (!opts.onlytheory)
      	datahistos[0].Draw((TH1F*)r_datatot->Clone(), "E3 same");
    }
  r_datatot->GetXaxis()->SetRange(datahistos[0].getlowrange(), datahistos[0].getuprange());

  //plot lines at 1 (or 0 for diff plots)
  TLine *r_ref;
  if (opts.diff)
    r_ref = new TLine(r_templ->GetBinLowEdge(r_templ->GetXaxis()->GetFirst()), 0, r_templ->GetXaxis()->GetBinUpEdge(r_templ->GetXaxis()->GetLast()), 0);
  else
    r_ref = new TLine(r_templ->GetBinLowEdge(r_templ->GetXaxis()->GetFirst()), 1, r_templ->GetXaxis()->GetBinUpEdge(r_templ->GetXaxis()->GetLast()), 1);
  r_ref->SetLineStyle(2);
  r_ref->SetLineStyle(1);
  r_ref->Draw();

  TLegend * legr;
  if (opts.multitheory)
    legr = new TLegend(lmarg+0.04, marg0/ry+0.03, lmarg+0.04+0.30, marg0/ry+0.03+0.04/ry);

  //Draw ratios
  for (vector <Subplot>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      (*it).getrthshift()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
      (*it).getrthshift()->SetLineStyle(2);
      (*it).getrthshift()->SetLineWidth(opts.lwidth);

      (*it).getrth()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
      (*it).getrth()->SetLineWidth(opts.lwidth);
      if (opts.bw)
	(*it).getrth()->SetLineStyle(opts.lstyles[labels[it-datahistos.begin()]]);

      vector <range> rthranges = historanges((*it).getrthshift());
      if (!opts.threepanels)
	{
	  for (vector<range>::iterator r = rthranges.begin(); r != rthranges.end(); r++)
	    {
	      (*it).getrthshift()->SetAxisRange((*r).lowedge, (*r).upedge);
	      if (!opts.onlytheory)
		if (!opts.multitheory || (it - datahistos.begin() == 0)) //if in multitheory mode, plot only the first theory
		  if (!opts.nothshifts)
		    if ((!opts.points || (*it).bincenter()) && (*it).nbins() > 1)
		      (*it).Draw((TH1F*)(*it).getrthshift()->Clone(), "LX same");
		    else
		      (*it).Draw((TH1F*)(*it).getrthshift()->Clone(), "hist ][ same");

	    }
	  (*it).getrthshift()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
	}
      
      if ((!opts.points || (*it).bincenter()) && (*it).nbins() > 1) //plot as continous line with dashed error bands
	{
	  for (vector<range>::iterator r = rthranges.begin(); r != rthranges.end(); r++)
	    {
	      (*it).getrth()->SetAxisRange((*r).lowedge, (*r).upedge);
	      if (!opts.multitheory || (it - datahistos.begin() == 0)) //if in multitheory mode, plot only the first theory
		(*it).Draw((TH1F*)(*it).getrth()->Clone(), "LX same");
	    }
	  (*it).getrth()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
	  if (opts.therr)
	    {
	      (*it).getrtherr()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	      (*it).getrtherr()->SetMarkerSize(0);
	      (*it).getrtherr()->SetFillColor(opts.colors[labels[it-datahistos.begin()]]);
	      (*it).getrtherr()->SetFillStyle(opts.styles[labels[it-datahistos.begin()]]);
	      float toterr = 0;
	      for (int b = 1; b <= (*it).gettherr()->GetNbinsX(); b++)
		toterr += (*it).gettherr()->GetBinError(b);
	      if (toterr > 0)
		{
		  for (vector<range>::iterator r = rthranges.begin(); r != rthranges.end(); r++)
		    {
		      (*it).getrtherr()->SetAxisRange((*r).lowedge, (*r).upedge);
		      if (!opts.multitheory || (it - datahistos.begin() == 0)) //if in multitheory mode, plot only the first theory
			(*it).Draw((TH1F*)(*it).getrtherr()->Clone(), "E3L same");
		    }
		  (*it).getrtherr()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
		}
	    }
	}
      else //plot as displaced TGraphs
	{
	  TGraphAsymmErrors * r_gtherr = new TGraphAsymmErrors((*it).getrth());
	  r_gtherr->SetMarkerStyle(opts.markers[labels[it-datahistos.begin()]]);
	  r_gtherr->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	  r_gtherr->SetMarkerSize(2 * opts.resolution / 1200);
	  r_gtherr->SetMarkerColor(opts.colors[labels[it-datahistos.begin()]]);
	  for (int b = 0; b < r_gtherr->GetN(); b++)
	    {
	      //Set X error to 0
	      r_gtherr->SetPointEXlow(b, 0);
	      r_gtherr->SetPointEXhigh(b, 0);

	      //displace horizontally
	      double x, y;
	      r_gtherr->GetPoint(b, x, y);
	      float width = (*it).getrth()->GetBinWidth(b + 1);
	      float lowedge = (*it).getrth()->GetBinLowEdge(b + 1);
	      x = lowedge + (it - datahistos.begin() + 1) * width/(datahistos.size() + 1);
	      r_gtherr->SetPoint(b, x, y);
	      //Set Y error
	      float errup, errdown;
	      if (opts.therr)
		{    
		  errup = (*it).getrtherrup()->GetBinContent(b + 1) - (*it).getrth()->GetBinContent(b + 1);
		  errdown = (*it).getrth()->GetBinContent(b + 1) - (*it).getrtherrdown()->GetBinContent(b + 1);
		}
	      else
		{
		  errup = 0;
		  errdown = 0;
		}
	      r_gtherr->SetPointEYhigh(b, errup);
	      r_gtherr->SetPointEYlow(b, errdown);
	    }
	  if (!opts.multitheory || (it - datahistos.begin() == 0)) //if in multitheory mode, plot only the first theory
	    r_gtherr->Draw("P same");
	}
      if (opts.multitheory && (it - datahistos.begin() == 0))
	{
	  if ((!opts.points || (*it).bincenter()) && (*it).nbins() > 1)
	    if (opts.therr && (*it).HasTherr())
	      legr->AddEntry((*it).getrtherr(), (labels[it-datahistos.begin()]).c_str(), "lf");
	    else
	      legr->AddEntry((*it).getrth(), (labels[it-datahistos.begin()]).c_str(), "l");
	  else
	    if (opts.therr && (*it).HasTherr())
	      legr->AddEntry((*it).getrth(), (labels[it-datahistos.begin()]).c_str(), "pe");
	    else
	      leg2->AddEntry((*it).getrth(), (labels[it-datahistos.begin()]).c_str(), "p");
	}
    }
  
  //draw theory error borders
  if (opts.therr)
    for (vector <Subplot>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
      {
	(*it).getrtherrup()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	(*it).getrtherrdown()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	(*it).getrtherrup()->SetLineWidth(opts.lwidth);
	(*it).getrtherrdown()->SetLineWidth(opts.lwidth);
	if ((!opts.points || (*it).bincenter()) && (*it).nbins() > 1)
	  {
	    vector <range> rthranges = historanges((*it).getth());
	    for (vector<range>::iterator r = rthranges.begin(); r != rthranges.end(); r++)
	      {
		(*it).getrtherrup()->SetAxisRange((*r).lowedge, (*r).upedge);
		if (!opts.multitheory || (it - datahistos.begin() == 0)) //if in multitheory mode, plot only the first theory
		  (*it).Draw((TH1F*)(*it).getrtherrup()->Clone(), "LX same");
		(*it).getrtherrdown()->SetAxisRange((*r).lowedge, (*r).upedge);
		if (!opts.multitheory || (it - datahistos.begin() == 0)) //if in multitheory mode, plot only the first theory
		  (*it).Draw((TH1F*)(*it).getrtherrdown()->Clone(), "LX same");
	      }
	    (*it).getrtherrup()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
	    (*it).getrtherrdown()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
	  }
      }

  //Draw data points
  if (!opts.onlytheory)
    datahistos[0].Draw(r_data, "PE1 same");

  if (opts.multitheory)
    {
      legr->SetFillColor(0);
      legr->SetBorderSize(0);
      legr->SetTextAlign(12);
      legr->SetTextFont(62);
      legr->SetTextSize(txtsize * 0.8/ry);
      legr->Draw();
    }

  //Theory+shifts/Data ratio Pad (optional)
  if (opts.multitheory && datahistos.size() == 3)
    {
      Shifts->cd();
      if (datahistos[0].getlogx())
	Shifts->SetLogx();

      //Set up template histogram for axis
      r_templ->GetYaxis()->SetLabelSize(txtsize/sy);
      r_templ->GetYaxis()->SetTitleSize(txtsize/sy);
      r_templ->GetYaxis()->SetTitleOffset((offset+0.3) * sy);

      //draw axis
      r_templ->DrawCopy("AXIS");

      //draw data
      vector <range> rdtranges = historanges(r_datatot);
      for (vector<range>::iterator r = rdtranges.begin(); r != rdtranges.end(); r++)
	{
	  r_datatot->SetAxisRange((*r).lowedge, (*r).upedge);
	  if (!opts.onlytheory)
	    datahistos[0].Draw((TH1F*)r_datatot->Clone(), "E3 same");
	}
      r_datatot->GetXaxis()->SetRange(datahistos[0].getlowrange(), datahistos[0].getuprange());

      //plot lines at 1 (or 0 for diff plots)
      TLine *r_ref;
      if (opts.diff)
	r_ref = new TLine(r_templ->GetBinLowEdge(r_templ->GetXaxis()->GetFirst()), 0, r_templ->GetXaxis()->GetBinUpEdge(r_templ->GetXaxis()->GetLast()), 0);
      else
	r_ref = new TLine(r_templ->GetBinLowEdge(r_templ->GetXaxis()->GetFirst()), 1, r_templ->GetXaxis()->GetBinUpEdge(r_templ->GetXaxis()->GetLast()), 1);
      r_ref->SetLineStyle(2);
      r_ref->SetLineStyle(1);
      r_ref->Draw();

      TLegend * legr = new TLegend(lmarg+0.04, marg0/sy+0.03, lmarg+0.04+0.30, marg0/sy+0.03+0.04/sy);

      //Draw ratios
      for (vector <Subplot>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
	{
	  (*it).getrthshift()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	  (*it).getrthshift()->SetLineStyle(2);
	  (*it).getrthshift()->SetLineWidth(opts.lwidth);

	  (*it).getrth()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	  (*it).getrth()->SetLineWidth(opts.lwidth);
	  if (opts.bw)
	    (*it).getrth()->SetLineStyle(opts.lstyles[labels[it-datahistos.begin()]]);

	  vector <range> rthranges = historanges((*it).getrthshift());
	  if (!opts.threepanels)
	    {
	      for (vector<range>::iterator r = rthranges.begin(); r != rthranges.end(); r++)
		{
		  (*it).getrthshift()->SetAxisRange((*r).lowedge, (*r).upedge);
		  if (!opts.onlytheory)
		    if (!opts.multitheory || (it - datahistos.begin() == 1)) //if in multitheory mode, plot only the second theory
		      if (!opts.nothshifts)
			if ((!opts.points || (*it).bincenter()) && (*it).nbins() > 1)
			  (*it).Draw((TH1F*)(*it).getrthshift()->Clone(), "LX same");
			else
			  (*it).Draw((TH1F*)(*it).getrthshift()->Clone(), "hist ][ same");
		}
	      (*it).getrthshift()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
	    }
      
	  if ((!opts.points || (*it).bincenter()) && (*it).nbins() > 1) //plot as continous line with dashed error bands
	    {
	      for (vector<range>::iterator r = rthranges.begin(); r != rthranges.end(); r++)
		{
		  (*it).getrth()->SetAxisRange((*r).lowedge, (*r).upedge);
		  if (!opts.multitheory || (it - datahistos.begin() == 1)) //if in multitheory mode, plot only the second theory
		    (*it).Draw((TH1F*)(*it).getrth()->Clone(), "LX same");
		}
	      (*it).getrth()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
	      if (opts.therr)
		{
		  (*it).getrtherr()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
		  (*it).getrtherr()->SetMarkerSize(0);
		  (*it).getrtherr()->SetFillColor(opts.colors[labels[it-datahistos.begin()]]);
		  (*it).getrtherr()->SetFillStyle(opts.styles[labels[it-datahistos.begin()]]);
		  float toterr = 0;
		  for (int b = 1; b <= (*it).gettherr()->GetNbinsX(); b++)
		    toterr += (*it).gettherr()->GetBinError(b);
		  if (toterr > 0)
		    {
		      for (vector<range>::iterator r = rthranges.begin(); r != rthranges.end(); r++)
			{
			  (*it).getrtherr()->SetAxisRange((*r).lowedge, (*r).upedge);
			  if (!opts.multitheory || (it - datahistos.begin() == 1)) //if in multitheory mode, plot only the second theory
			    (*it).Draw((TH1F*)(*it).getrtherr()->Clone(), "E3L same");
			}
		      (*it).getrtherr()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
		    }
		}
	    }
	  else //plot as displaced TGraphs
	    {
	      TGraphAsymmErrors * r_gtherr = new TGraphAsymmErrors((*it).getrth());
	      r_gtherr->SetMarkerStyle(opts.markers[labels[it-datahistos.begin()]]);
	      r_gtherr->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	      r_gtherr->SetMarkerSize(2 * opts.resolution / 1200);
	      r_gtherr->SetMarkerColor(opts.colors[labels[it-datahistos.begin()]]);
	      for (int b = 0; b < r_gtherr->GetN(); b++)
		{
		  //Set X error to 0
		  r_gtherr->SetPointEXlow(b, 0);
		  r_gtherr->SetPointEXhigh(b, 0);

		  //displace horizontally
		  double x, y;
		  r_gtherr->GetPoint(b, x, y);
		  float width = (*it).getrth()->GetBinWidth(b + 1);
		  float lowedge = (*it).getrth()->GetBinLowEdge(b + 1);
		  x = lowedge + (it - datahistos.begin() + 1) * width/(datahistos.size() + 1);
		  r_gtherr->SetPoint(b, x, y);
		  //Set Y error
		  float errup, errdown;
		  if (opts.therr)
		    {    
		      errup = (*it).getrtherrup()->GetBinContent(b + 1) - (*it).getrth()->GetBinContent(b + 1);
		      errdown = (*it).getrth()->GetBinContent(b + 1) - (*it).getrtherrdown()->GetBinContent(b + 1);
		    }
		  else
		    {
		      errup = 0;
		      errdown = 0;
		    }
		  r_gtherr->SetPointEYhigh(b, errup);
		  r_gtherr->SetPointEYlow(b, errdown);
		}
	      if (!opts.multitheory || (it - datahistos.begin() == 1)) //if in multitheory mode, plot only the second theory
		r_gtherr->Draw("P same");
	    }
	  if (it - datahistos.begin() == 1)
	    {
	      if ((!opts.points || (*it).bincenter()) && (*it).nbins() > 1) 
		if (opts.therr && (*it).HasTherr())
		  legr->AddEntry((*it).getrtherr(), (labels[it-datahistos.begin()]).c_str(), "lf");
		else
		  legr->AddEntry((*it).getrth(), (labels[it-datahistos.begin()]).c_str(), "l");
	      else
		if (opts.therr && (*it).HasTherr())
		  legr->AddEntry((*it).getrth(), (labels[it-datahistos.begin()]).c_str(), "pe");
		else
		  leg2->AddEntry((*it).getrth(), (labels[it-datahistos.begin()]).c_str(), "p");
	    }
	}
  
      //draw theory error borders
      if (opts.therr)
	for (vector <Subplot>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
	  {
	    (*it).getrtherrup()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	    (*it).getrtherrdown()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	    (*it).getrtherrup()->SetLineWidth(opts.lwidth);
	    (*it).getrtherrdown()->SetLineWidth(opts.lwidth);
	    if ((!opts.points || (*it).bincenter()) && (*it).nbins() > 1) 
	      {
		vector <range> rthranges = historanges((*it).getth());
		for (vector<range>::iterator r = rthranges.begin(); r != rthranges.end(); r++)
		  {
		    (*it).getrtherrup()->SetAxisRange((*r).lowedge, (*r).upedge);
		    if (!opts.multitheory || (it - datahistos.begin() == 1)) //if in multitheory mode, plot only the second theory
		      (*it).Draw((TH1F*)(*it).getrtherrup()->Clone(), "LX same");
		    (*it).getrtherrdown()->SetAxisRange((*r).lowedge, (*r).upedge);
		    if (!opts.multitheory || (it - datahistos.begin() == 1)) //if in multitheory mode, plot only the second theory
		      (*it).Draw((TH1F*)(*it).getrtherrdown()->Clone(), "LX same");
		  }
		(*it).getrtherrup()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
		(*it).getrtherrdown()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
	      }
	  }

      //Draw data points
      if (!opts.onlytheory)
	datahistos[0].Draw(r_data, "PE1 same");

      legr->SetFillColor(0);
      legr->SetBorderSize(0);
      legr->SetTextAlign(12);
      legr->SetTextFont(62);
      legr->SetTextSize(txtsize * 0.8/sy);
      legr->Draw();
    }
  else if (opts.threepanels)
    {
      Shifts->cd();
      if (datahistos[0].getlogx())
	Shifts->SetLogx();

      //Set up template histogram for axis
      r_templ->GetYaxis()->SetLabelSize(txtsize/sy);
      r_templ->GetYaxis()->SetTitleSize(txtsize/sy);
      r_templ->GetYaxis()->SetTitleOffset((offset+0.3) * sy);
      string ytitle = "";
      if (opts.diff)
	{
	  if (opts.ratiototheory)
	    ytitle = (string) "Data-" + opts.theorylabel;
	  else
	    ytitle = "Theory-Data";
	}
      else
	{
	  if (opts.ratiototheory)
	    ytitle = (string) "Ratio to " + opts.theorylabel;
	  else
	    ytitle = "#frac{Theory+shifts}{Data}";
	}

      r_templ->SetYTitle(ytitle.c_str());

      //Evaluate maximum and minimum
      mx = 0;
      TH1F * r_dataerr = (TH1F*) r_data->Clone();
      for (int b = 1; b <= r_data->GetNbinsX(); b++)
	r_dataerr->SetBinContent(b, r_data->GetBinContent(b) + r_data->GetBinError(b));
      if (!opts.onlytheory)
	mx = r_dataerr->GetBinContent(r_dataerr->GetMaximumBin());
      for (vector <Subplot>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
	mx = max(mx, (float)((*it).getrthshift()->GetMaximum()));

      mn = mx;
      for (int b = 1; b <= r_dataerr->GetNbinsX(); b++)
	r_dataerr->SetBinContent(b, r_data->GetBinContent(b) - r_data->GetBinError(b));
      if (!opts.onlytheory)
	mn = hmin(r_dataerr);
      for (vector <Subplot>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
	mn = min(mn, (float)(hmin((*it).getrthshift())));
      float delta = mx - mn;
      if (datahistos[0].getymaxr() != 0)
	{
	  mx = datahistos[0].getymaxr();
	  mn = datahistos[0].getyminr();
	  delta = 0;
	}
      r_templ->SetMaximum(mx + delta * 0.2);
      r_templ->SetMinimum(mn - delta * 0.2);

      r_templ->DrawCopy("AXIS");

      /*
      //plot data
      if (!opts.onlytheory)
	datahistos[0].Draw(r_data, "PE1 same");
      */

      //plot lines at 1 (or 0 for diff plots)
      TLine *rs_ref;
      if (opts.diff)
	rs_ref = new TLine(r_templ->GetBinLowEdge(r_templ->GetXaxis()->GetFirst()), 0, r_templ->GetXaxis()->GetBinUpEdge(r_templ->GetXaxis()->GetLast()), 0);
      else
	rs_ref = new TLine(r_templ->GetBinLowEdge(r_templ->GetXaxis()->GetFirst()), 1, r_templ->GetXaxis()->GetBinUpEdge(r_templ->GetXaxis()->GetLast()), 1);
      rs_ref->SetLineStyle(2);
      rs_ref->SetLineStyle(1);
      rs_ref->Draw();

      //Draw ratios
      for (vector <Subplot>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
	{
	  vector <range> rthranges = historanges((*it).getrthshift());
	  for (vector<range>::iterator r = rthranges.begin(); r != rthranges.end(); r++)
	    {
	      (*it).getrthshift()->SetAxisRange((*r).lowedge, (*r).upedge);
	      if (!opts.nothshifts)
		if ((!opts.points || (*it).bincenter()) && (*it).nbins() > 1)
		  (*it).Draw((TH1F*)(*it).getrthshift()->Clone(), "LX same");
		else
		  (*it).Draw((TH1F*)(*it).getrthshift()->Clone(), "hist ][ same");
	    }
	  (*it).getrthshift()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
	}	  
      //plot data
      if (!opts.onlytheory)
	datahistos[0].Draw(r_data, "PE1 same");
    }

  //Theory-Data pulls pad
  if (opts.multitheory && datahistos.size() >= 2)
    {
      int nth = datahistos.size();
      Pulls->cd();
      if (datahistos[0].getlogx())
	Pulls->SetLogx();

      //Set up template histogram for axis
      r_templ->GetXaxis()->SetLabelSize(txtsize/py);
      r_templ->GetXaxis()->SetTitleSize(txtsize/py);

      r_templ->GetYaxis()->SetLabelSize(txtsize/py);
      r_templ->GetYaxis()->SetTitleSize(txtsize/py);
      r_templ->GetYaxis()->SetTitleOffset((offset+0.3) * py);

      //draw axis
      r_templ->DrawCopy("AXIS");

      //draw data
      vector <range> rdtranges = historanges(r_datatot);
      for (vector<range>::iterator r = rdtranges.begin(); r != rdtranges.end(); r++)
	{
	  r_datatot->SetAxisRange((*r).lowedge, (*r).upedge);
	  if (!opts.onlytheory)
	    datahistos[0].Draw((TH1F*)r_datatot->Clone(), "E3 same");
	}
      r_datatot->GetXaxis()->SetRange(datahistos[0].getlowrange(), datahistos[0].getuprange());

      //plot lines at 1 (or 0 for diff plots)
      TLine *r_ref;
      if (opts.diff)
	r_ref = new TLine(r_templ->GetBinLowEdge(r_templ->GetXaxis()->GetFirst()), 0, r_templ->GetXaxis()->GetBinUpEdge(r_templ->GetXaxis()->GetLast()), 0);
      else
	r_ref = new TLine(r_templ->GetBinLowEdge(r_templ->GetXaxis()->GetFirst()), 1, r_templ->GetXaxis()->GetBinUpEdge(r_templ->GetXaxis()->GetLast()), 1);
      r_ref->SetLineStyle(2);
      r_ref->SetLineStyle(1);
      r_ref->Draw();

      TLegend * legr = new TLegend(lmarg+0.04, bmarg/py+0.03, lmarg+0.04+0.30, bmarg/py+0.03+0.04/py);

      //Draw ratios
      for (vector <Subplot>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
	{
	  (*it).getrthshift()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	  (*it).getrthshift()->SetLineStyle(2);
	  (*it).getrthshift()->SetLineWidth(opts.lwidth);

	  (*it).getrth()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	  (*it).getrth()->SetLineWidth(opts.lwidth);
	  if (opts.bw)
	    (*it).getrth()->SetLineStyle(opts.lstyles[labels[it-datahistos.begin()]]);

	  vector <range> rthranges = historanges((*it).getrthshift());
	  if (!opts.threepanels)
	    {
	      for (vector<range>::iterator r = rthranges.begin(); r != rthranges.end(); r++)
		{
		  (*it).getrthshift()->SetAxisRange((*r).lowedge, (*r).upedge);
		  if (!opts.onlytheory)
		    if (!opts.multitheory || (it - datahistos.begin() == nth-1)) //if in multitheory mode, plot only the second theory
		      if (!opts.nothshifts)
			if ((!opts.points || (*it).bincenter()) && (*it).nbins() > 1) 
			  (*it).Draw((TH1F*)(*it).getrthshift()->Clone(), "LX same");
			else
			  (*it).Draw((TH1F*)(*it).getrthshift()->Clone(), "hist ][ same");
		}
	      (*it).getrthshift()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
	    }
      
	  if ((!opts.points || (*it).bincenter()) && (*it).nbins() > 1) //plot as continous line with dashed error bands
	    {
	      for (vector<range>::iterator r = rthranges.begin(); r != rthranges.end(); r++)
		{
		  (*it).getrth()->SetAxisRange((*r).lowedge, (*r).upedge);
		  if (!opts.multitheory || (it - datahistos.begin() == nth-1)) //if in multitheory mode, plot only the second theory
		    (*it).Draw((TH1F*)(*it).getrth()->Clone(), "LX same");
		}
	      (*it).getrth()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
	      if (opts.therr)
		{
		  (*it).getrtherr()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
		  (*it).getrtherr()->SetMarkerSize(0);
		  (*it).getrtherr()->SetFillColor(opts.colors[labels[it-datahistos.begin()]]);
		  (*it).getrtherr()->SetFillStyle(opts.styles[labels[it-datahistos.begin()]]);
		  float toterr = 0;
		  for (int b = 1; b <= (*it).gettherr()->GetNbinsX(); b++)
		    toterr += (*it).gettherr()->GetBinError(b);
		  if (toterr > 0)
		    {
		      for (vector<range>::iterator r = rthranges.begin(); r != rthranges.end(); r++)
			{
			  (*it).getrtherr()->SetAxisRange((*r).lowedge, (*r).upedge);
			  if (!opts.multitheory || (it - datahistos.begin() == nth-1)) //if in multitheory mode, plot only the second theory
			    (*it).Draw((TH1F*)(*it).getrtherr()->Clone(), "E3L same");
			}
		      (*it).getrtherr()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
		    }
		}
	    }
	  else //plot as displaced TGraphs
	    {
	      TGraphAsymmErrors * r_gtherr = new TGraphAsymmErrors((*it).getrth());
	      r_gtherr->SetMarkerStyle(opts.markers[labels[it-datahistos.begin()]]);
	      r_gtherr->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	      r_gtherr->SetMarkerSize(2 * opts.resolution / 1200);
	      r_gtherr->SetMarkerColor(opts.colors[labels[it-datahistos.begin()]]);
	      for (int b = 0; b < r_gtherr->GetN(); b++)
		{
		  //Set X error to 0
		  r_gtherr->SetPointEXlow(b, 0);
		  r_gtherr->SetPointEXhigh(b, 0);

		  //displace horizontally
		  double x, y;
		  r_gtherr->GetPoint(b, x, y);
		  float width = (*it).getrth()->GetBinWidth(b + 1);
		  float lowedge = (*it).getrth()->GetBinLowEdge(b + 1);
		  x = lowedge + (it - datahistos.begin() + 1) * width/(datahistos.size() + 1);
		  r_gtherr->SetPoint(b, x, y);
		  //Set Y error
		  float errup, errdown;
		  if (opts.therr)
		    {    
		      errup = (*it).getrtherrup()->GetBinContent(b + 1) - (*it).getrth()->GetBinContent(b + 1);
		      errdown = (*it).getrth()->GetBinContent(b + 1) - (*it).getrtherrdown()->GetBinContent(b + 1);
		    }
		  else
		    {
		      errup = 0;
		      errdown = 0;
		    }
		  r_gtherr->SetPointEYhigh(b, errup);
		  r_gtherr->SetPointEYlow(b, errdown);
		}
	      if (!opts.multitheory || (it - datahistos.begin() == nth-1)) //if in multitheory mode, plot only the second theory
		r_gtherr->Draw("P same");
	    }
	  if (it - datahistos.begin() == nth-1)
	    {
	      if ((!opts.points || (*it).bincenter()) && (*it).nbins() > 1) 
		if (opts.therr && (*it).HasTherr())
		  legr->AddEntry((*it).getrtherr(), (labels[it-datahistos.begin()]).c_str(), "lf");
		else
		  legr->AddEntry((*it).getrth(), (labels[it-datahistos.begin()]).c_str(), "l");
	      else
		if (opts.therr && (*it).HasTherr())
		  legr->AddEntry((*it).getrth(), (labels[it-datahistos.begin()]).c_str(), "pe");
		else
		  leg2->AddEntry((*it).getrth(), (labels[it-datahistos.begin()]).c_str(), "p");
	    }
	}	  
  
      //draw theory error borders
      if (opts.therr)
	for (vector <Subplot>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
	  {
	    (*it).getrtherrup()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	    (*it).getrtherrdown()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	    (*it).getrtherrup()->SetLineWidth(opts.lwidth);
	    (*it).getrtherrdown()->SetLineWidth(opts.lwidth);
	    if ((!opts.points || (*it).bincenter()) && (*it).nbins() > 1) 
	      {
		vector <range> rthranges = historanges((*it).getth());
		for (vector<range>::iterator r = rthranges.begin(); r != rthranges.end(); r++)
		  {
		    (*it).getrtherrup()->SetAxisRange((*r).lowedge, (*r).upedge);
		    if (!opts.multitheory || (it - datahistos.begin() == nth-1)) //if in multitheory mode, plot only the second theory
		      (*it).Draw((TH1F*)(*it).getrtherrup()->Clone(), "LX same");
		    (*it).getrtherrdown()->SetAxisRange((*r).lowedge, (*r).upedge);
		    if (!opts.multitheory || (it - datahistos.begin() == nth-1)) //if in multitheory mode, plot only the second theory
		      (*it).Draw((TH1F*)(*it).getrtherrdown()->Clone(), "LX same");
		  }
		(*it).getrtherrup()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
		(*it).getrtherrdown()->GetXaxis()->SetRange((*it).getlowrange(), (*it).getuprange());
	      }
	  }

      //Draw data points
      if (!opts.onlytheory)
	datahistos[0].Draw(r_data, "PE1 same");

      legr->SetFillColor(0);
      legr->SetBorderSize(0);
      legr->SetTextAlign(12);
      legr->SetTextFont(62);
      legr->SetTextSize(txtsize * 0.8/py);
      legr->Draw();
    }
  else if (opts.twopanels || opts.threepanels)
    {
      Pulls->cd();
      if (datahistos[0].getlogx())
	Pulls->SetLogx();

      TH1F * pull = datahistos[0].getpull();

      r_templ->GetYaxis()->SetLabelSize(txtsize/py);
      r_templ->GetYaxis()->SetTitleSize(txtsize/py);
      r_templ->GetYaxis()->SetTitleOffset((offset+0.3) * py);

      //r_templ->SetYTitle("#frac{Theory+shifts - Data}{#sigma uncor}");
      r_templ->SetYTitle("pulls   ");

      r_templ->GetXaxis()->SetLabelSize(txtsize/py);
      r_templ->GetXaxis()->SetTitleSize(txtsize/py);
      //  r_templ->GetXaxis()->SetTitleOffset(1);

      r_templ->SetStats(0);
      r_templ->SetMinimum(-3.5);
      r_templ->SetMaximum(3.5);

      //plot axis
      r_templ->Draw("AXIS");

      //plot lines at 1, -1, 0
      TLine *one = new TLine(r_templ->GetBinLowEdge(r_templ->GetXaxis()->GetFirst()), 1, r_templ->GetXaxis()->GetBinUpEdge(r_templ->GetXaxis()->GetLast()), 1);
      one->SetLineStyle(2);
      TLine *minusone = new TLine(r_templ->GetBinLowEdge(r_templ->GetXaxis()->GetFirst()), -1, r_templ->GetXaxis()->GetBinUpEdge(r_templ->GetXaxis()->GetLast()), -1);
      minusone->SetLineStyle(2);
      TLine *zero = new TLine(r_templ->GetBinLowEdge(r_templ->GetXaxis()->GetFirst()), 0, r_templ->GetXaxis()->GetBinUpEdge(r_templ->GetXaxis()->GetLast()), 0);
      one->Draw();
      minusone->Draw();
      zero->Draw();

      //plot pulls
      for (vector <Subplot>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
	{
	  if (datahistos.size() == 1)
	    {
	      (*it).getpull()->SetFillColor(opts.colors[labels[it-datahistos.begin()]]);
	      (*it).getpull()->SetFillStyle(1001);
	    }
	  (*it).getpull()->SetLineStyle(1);
	  (*it).getpull()->SetLineWidth(opts.lwidth);
	  (*it).getpull()->SetLineColor(opts.colors[labels[it-datahistos.begin()]]);
	  datahistos[0].Draw((TH1F*)(*it).getpull()->Clone(), "same ][");
	}	  
    }

  //Labels
  if (opts.twopanels || opts.threepanels)
    {
      cnv->cd(1);
      DrawLabels();
      if (opts.drawlogo)
        DrawLogo()->Draw();   
      cnv->cd(2);
      DrawLabels();
    }
  else
    {
      cnv->cd();
      DrawLabels();
      if (opts.drawlogo)
        DrawLogo()->Draw();   
    }

  return cnv;
}
