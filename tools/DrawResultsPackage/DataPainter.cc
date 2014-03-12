#include <DataPainter.h>
#include <CommandParser.h>

#include <DrawLogo.h>

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
dataseth::dataseth(string dataname, string dir, string lab, DataSet* DT, int subp, int colindx): name(dataname), label(lab), coli(colindx)
{
  TGraphErrors * datagraph =          DT->GetDataTot(subp);	       
  vector <float> bins1 =              DT->getbins1(subp);	       
  vector <float> bins2 =	      DT->getbins2(subp);	       
  vector <float> data = 	      DT->getdata(subp);	       
  vector <float> uncorerr = 	      DT->getuncor(subp);	       
  vector <float> toterr = 	      DT->gettoterr(subp);	       
  vector <float> theory = 	      DT->gettheory(subp);	       
  vector <float> theoryshifted =      DT->gettheoryshifted(subp);   
  vector <float> therrup = 	      DT->gettherrup(subp);	       
  vector <float> therrdown = 	      DT->gettherrdown(subp);       
  vector <float> pulls = 	      DT->getpulls(subp);	       
  logx = 			      DT->GetXlog(subp);	       
  logy = 			      DT->GetYlog(subp);	       
  xmin = 			      DT->GetXmin(subp);	       
  xmax =			      DT->GetXmax(subp);	       
  yminr = 			      DT->GetYminR(subp);	       
  ymaxr =			      DT->GetYmaxR(subp);	       
  string xlabel = 		      DT->GetXTitle(subp);	       
  string ylabel =		      DT->GetYTitle(subp);         
  title =		              DT->GetTitle(subp);
  extralabel =		              DT->GetLabel(subp);
  experiment =		              DT->GetExperiment(subp);

  //bins array
  float bin[bins1.size() + 1];
  vector <double> bins;

  //If no bin edges information is available, need to use TGraph for plotting
  maketgraph = false;
  for (int i = 0; i < datagraph->GetN(); i++)
    if(datagraph->GetX()[i] != 0)
      maketgraph = true;

  if (maketgraph)
    {
      //save original x values
      for (int i = 0; i < datagraph->GetN(); i++)
	valx.push_back(datagraph->GetX()[i]);

      //make arbitrary bin edges
      vector <double> temp;
      for (int i = 0; i < datagraph->GetN(); i++)
	temp.push_back(datagraph->GetX()[i]);

      sort(temp.begin(), temp.end());

      //adjust x axis
      double xaxmin = *(temp.begin());
      double xaxmax = *(temp.end()-1);
      if (datagraph->GetN() > 1 && (*(temp.end()-1) - *(temp.begin()) > 0))
	{
	  if (logx)
	    {
	      double axislength = *(temp.end()-1) / *(temp.begin());
	      xaxmin = xaxmin / pow(axislength,1./20);
	      xaxmax = xaxmax * pow(axislength,1./20);
	    }
	  else
	    {
	      double axislength = *(temp.end()-1) - *(temp.begin());
	      xaxmin = xaxmin - axislength/20.;
	      xaxmax = xaxmax + axislength/20.;
	    }
	}
      else
	{
	  xaxmin *= 0.9999;
	  xaxmax *= 1.0001;
	}
      bins.push_back(xaxmin);
      if (datagraph->GetN() > 1)
	for (vector<double>::iterator it = temp.begin()+1; it != temp.end(); it++)
	  bins.push_back((*it + *(it-1))/2);
      bins.push_back(xaxmax);


      for (vector<double>::iterator it = bins.begin(); it != bins.end(); it++)
	bin[it-bins.begin()] = *it;
    }
  else
    {
      //fill empty gap
      int pos = 0;
      float bmin, bmax;
      while (pos != -1)
	{
	  pos = -1;
	  vector<float>::iterator it1 = bins1.begin();
	  vector<float>::iterator it2 = bins2.begin();
	  for (; (it1+1) != bins1.end(); it1++, it2++)
	    if (*(it1+1) != *it2 && *it2 < *(it1+1))
	      {
		pos = (it1 - bins1.begin()) + 1;
		bmin = *it2;
		bmax = *(it1+1);
	      }
	  if (pos != -1)
	    {
	      bins1.insert(bins1.begin()+pos, bmin);
	      bins2.insert(bins2.begin()+pos, bmax);
	      data.insert(data.begin()+pos, 0);
	      uncorerr.insert(uncorerr.begin()+pos, 0); 
	      toterr.insert(toterr.begin()+pos, 0); 
	      theory.insert(theory.begin()+pos, 0); 
	      theoryshifted.insert(theoryshifted.begin() +pos, 0); 
	      therrup.insert(therrup.begin()+pos, 0); 
	      therrdown.insert(therrdown.begin()+pos, 0); 
	      pulls.insert(pulls.begin()+pos, 0);
	    }
	}

      //make bins array
      int i = 0;
      for (vector<float>::iterator it = bins1.begin(); it != bins1.end(); it++)
	{
	  bin[i] = *it;
	  i++;
	}
      bin[i] = *(bins2.end()-1);
    }

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
  else
    {
      xmin = hdata->GetXaxis()->GetBinLowEdge(hdata->GetXaxis()->GetFirst());
      xmax = hdata->GetXaxis()->GetBinUpEdge(hdata->GetXaxis()->GetLast() - 1);
    }
      
  hdata->SetXTitle(xlabel.c_str());
  hdata->SetYTitle(ylabel.c_str());

  //set rapidity as default label 
  if (xlabel == "" && ylabel == "")
    {
      hdata->SetXTitle("(Set XTitle:<label>)");
      hpull->SetXTitle("(Set XTitle:<label>)");
      hdata->SetYTitle("(Set YTitle:<label>)");
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
      if (!opts.ratiototheory)
	hpull->SetBinContent(b + 1, -pulls[b]);
      else
	hpull->SetBinContent(b + 1, pulls[b]);
      hpull->SetBinError(b + 1, 0);
    }

  //Prepare ratio histograms
  r_th = (TH1F*)hth->Clone();
  r_thshift = (TH1F*)hthshift->Clone();
  r_therr = (TH1F*)htherr->Clone();
  r_therrup = (TH1F*)htherrup->Clone();
  r_therrdown = (TH1F*)htherrdown->Clone();

  //Select reference, data or theory
  //  TH1F * refdata;
  if (opts.ratiototheory)
    href = (TH1F*)hth->Clone();
  else
    href = (TH1F*)hdata->Clone();

  //Prevent double counting of errors in ratio
  for (int b = 1; b <= href->GetNbinsX(); b++)
    href->SetBinError(b, 0);

  /*
  r_th->Divide(refdata);
  r_thshift->Divide(refdata);
  r_therr->Divide(refdata);
  r_therrup->Divide(refdata);
  r_therrdown->Divide(refdata);

  for (int b = 1; b <= r_th->GetNbinsX(); b++)
    r_th->SetBinError(b, 0);
  for (int b = 1; b <= r_thshift->GetNbinsX(); b++)
    r_thshift->SetBinError(b, 0);
  for (int b = 1; b <= r_therr->GetNbinsX(); b++)
    r_therr->SetBinError(b, (r_therrup->GetBinContent(b) - r_therrdown->GetBinContent(b)) / 2 );
  for (int b = 1; b <= r_therrup->GetNbinsX(); b++)
    r_therrup->SetBinError(b, 0);
  for (int b = 1; b <= r_therrdown->GetNbinsX(); b++)
    r_therrdown->SetBinError(b, 0);
  */
}
void dataseth::Draw(TH1F* histo, string opt)
{
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
	opt.erase(opt.find("E1") + 1);

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

TCanvas * DataPainter(int dataindex, vector <dataseth> datahistos)
{
  if (datahistos.size() < 1)
    {
      cout << "Empty dataset vector for dataset idx: " << dataindex << endl;
      return 0;
    }

  char cnvname[15];
  sprintf(cnvname, "%d_pulls",  dataindex);

  TCanvas * cnv;
  if (opts.twopanels || opts.threepanels)
    cnv = new TCanvas(cnvname, "", 0, 0, 2 * opts.resolution, opts.resolution);
  else
    cnv = new TCanvas(cnvname, "", 0, 0, opts.resolution, opts.resolution);
  cnv->cd();

  TH1F * data = datahistos[0].getdata();
  TH1F * datatot = datahistos[0].getdatatot();
  string dataname = datahistos[0].getname();

  //Set the pads geometry
  //panels height
  float dy;

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
  float my, ry, sy, py;
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
  TAxis * XAx = data->GetXaxis();
  float axmin = XAx->GetBinLowEdge(XAx->GetFirst());
  float axmax = XAx->GetBinUpEdge(XAx->GetLast());
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
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      if (opts.therr)
	mx = max(mx, (float)((*it).gettherrup()->GetMaximum()));
      else
	mx = max(mx, (float)((*it).getth()->GetMaximum()));
      mx = max(mx, (float)((*it).getthshift()->GetMaximum()));
    }

  float mn = mx;
  for (int b = 1; b <= datatot->GetNbinsX(); b++)
    dataerr->SetBinContent(b, datatot->GetBinContent(b) - datatot->GetBinError(b));
  if (!opts.onlytheory)
    mn = hmin(dataerr);
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      if (opts.therr)
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

  up_templ->SetMaximum(mx);
  up_templ->SetMinimum(mn);
  up_templ->SetStats(0);
  up_templ->Draw("AXIS");

  data->SetStats(0);
  data->SetLineColor(1);
  data->SetMarkerStyle(20);
  data->SetMarkerSize(2 * opts.resolution / 1200);
  if (!opts.onlytheory)
    datahistos[0].Draw(data, "PE1 same");

  datatot->SetFillColor(kYellow);
  datatot->SetLineColor(kYellow);
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
  datatot->SetAxisRange(datahistos[0].getxmin(), datahistos[0].getxmax());

  if (datahistos[0].getextralabel() != "")
    {
      TLatex l;
      l.SetNDC();
      l.SetTextFont(42);
      l.SetTextSize(0.04/my);
      l.DrawLatex(lmarg+0.05, (1-tmarg/my) - 0.13/my, datahistos[0].getextralabel().c_str());
    }

  //Main legend
  TLegend * leg = new TLegend(lmarg+0.02+0.02, mb+0.03, mb+0.4/my, lmarg+0.20);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextAlign(12);
  leg->SetTextSize(txtsize * 0.8/my);
  leg->SetTextFont(62);
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
  mark->SetMarkerStyle(opts.markers[0]);
  mark->SetMarkerSize(2 * opts.resolution / 1200);
  mark->SetMarkerColor(kBlack);
  TLine *cont = new TLine(0, 1, 1, 1);
  cont->SetLineStyle(1);
  TLine *dash = new TLine(0, 1, 1, 1);
  dash->SetLineStyle(2);
  if (opts.points && !datahistos[0].bincenter())
    leg->AddEntry(mark, opts.theorylabel.c_str(), "p");
  else
    leg->AddEntry(cont, opts.theorylabel.c_str(), "l");
  if (!opts.onlytheory)
    leg->AddEntry(dash, (opts.theorylabel + " + shifts").c_str(), "l");

  //Auxiliary legend
  TLegend * leg2 = new TLegend(lmarg+0.4, mb+0.03, 1-rmarg-0.01, mb+0.03 + datahistos.size() * 0.045/my);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextAlign(12);
  leg2->SetTextFont(62);
  leg2->SetTextSize(txtsize * 0.8/my);

  //Plot theories
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      TGraphAsymmErrors * gtherr = new TGraphAsymmErrors((*it).getth());
      (*it).getthshift()->SetLineColor(opts.colors[(*it).getcoli()]);
      (*it).getthshift()->SetLineStyle(2);
      (*it).getthshift()->SetLineWidth(2);

      vector <range> thranges = historanges((*it).getthshift());
      for (vector<range>::iterator r = thranges.begin(); r != thranges.end(); r++)
	{
	  (*it).getthshift()->SetAxisRange((*r).lowedge, (*r).upedge);
	  if (!opts.onlytheory)
	    (*it).Draw((TH1F*)(*it).getthshift()->Clone(), "LX same");
	}
      (*it).getthshift()->SetAxisRange((*it).getxmin(), (*it).getxmax());


      (*it).getth()->SetLineColor(opts.colors[(*it).getcoli()]);
      (*it).getth()->SetLineWidth(2);
      if (!opts.points || (*it).bincenter()) //plot as continous line with dashed error bands
	{
	  vector <range> thranges = historanges((*it).getth());
	  for (vector<range>::iterator r = thranges.begin(); r != thranges.end(); r++)
	    {
	      (*it).getth()->SetAxisRange((*r).lowedge, (*r).upedge);
	      (*it).Draw((TH1F*)(*it).getth()->Clone(), "LX same");
	    }
	  (*it).getth()->SetAxisRange((*it).getxmin(), (*it).getxmax());

	  if (opts.therr)
	    {    
	      (*it).gettherr()->SetLineColor(opts.colors[(*it).getcoli()]);
	      (*it).gettherr()->SetMarkerSize(0);
	      (*it).gettherr()->SetFillColor(opts.colors[(*it).getcoli()]);
	      (*it).gettherr()->SetFillStyle(opts.styles[(*it).getcoli()]);
	      float toterr = 0;
	      for (int b = 1; b <= (*it).gettherr()->GetNbinsX(); b++)
		toterr += (*it).gettherr()->GetBinError(b);
	      if (toterr > 0)
		{
		  vector <range> thranges = historanges((*it).getth());
		  for (vector<range>::iterator r = thranges.begin(); r != thranges.end(); r++)
		    {
		      (*it).gettherr()->SetAxisRange((*r).lowedge, (*r).upedge);
		      (*it).Draw((TH1F*)(*it).gettherr(), "E3L same");
		    }
		  (*it).gettherr()->SetAxisRange((*it).getxmin(), (*it).getxmax());
		}
	    }
	}
      else //plot as displaced points with vertical error line
	{
	  gtherr->SetMarkerStyle(opts.markers[(*it).getcoli()]);
	  gtherr->SetLineColor(opts.colors[(*it).getcoli()]);
	  gtherr->SetMarkerSize(2 * opts.resolution / 1200);
	  gtherr->SetMarkerColor(opts.colors[(*it).getcoli()]);
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
      if (!opts.points || (*it).bincenter())
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
  if (opts.therr)
    for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
      {
	(*it).gettherrup()->SetLineColor(opts.colors[(*it).getcoli()]);
	(*it).gettherrdown()->SetLineColor(opts.colors[(*it).getcoli()]);
	(*it).gettherrup()->SetLineWidth(2);
	(*it).gettherrdown()->SetLineWidth(2);
	if (!opts.points || (*it).bincenter())
	  {
	    vector <range> thranges = historanges((*it).getth());
	    for (vector<range>::iterator r = thranges.begin(); r != thranges.end(); r++)
	      {
		(*it).gettherrup()->SetAxisRange((*r).lowedge, (*r).upedge);
		(*it).Draw((TH1F*)(*it).gettherrup()->Clone(), "LX same");
		(*it).gettherrdown()->SetAxisRange((*r).lowedge, (*r).upedge);
		(*it).Draw((TH1F*)(*it).gettherrdown()->Clone(), "LX same");
	      }
	    (*it).gettherrup()->SetAxisRange((*it).getxmin(), (*it).getxmax());
	    (*it).gettherrdown()->SetAxisRange((*it).getxmin(), (*it).getxmax());
	  }
      }
  leg->Draw();
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
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
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
	  (*it).getrth()->Divide(refdata);
	  (*it).getrthshift()->Divide(refdata);
	  (*it).getrtherr()->Divide(refdata);
	  (*it).getrtherrup()->Divide(refdata);
	  (*it).getrtherrdown()->Divide(refdata);
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
      if (opts.ratiototheory)
	ytitle = (string) "Data-" + opts.theorylabel;
      else
	ytitle = "Theory-Data";
    }
  else
    {
      if (opts.ratiototheory)
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
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      if (opts.therr)
	mx = max(mx, (float)((*it).getrtherrup()->GetMaximum()));
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
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      if (opts.therr)
	mn = min(mn, (float)(hmin((*it).getrtherrdown())));
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
  r_datatot->SetAxisRange(datahistos[0].getxmin(), datahistos[0].getxmax());

  if (!opts.onlytheory)
    datahistos[0].Draw(r_data, "PE1 same");

  //plot lines at 1
  TLine *r_one = new TLine(r_templ->GetBinLowEdge(r_templ->GetXaxis()->GetFirst()), 1, r_templ->GetXaxis()->GetBinUpEdge(r_templ->GetXaxis()->GetLast()), 1);
  r_one->SetLineStyle(2);
  r_one->SetLineStyle(1);
  r_one->Draw();

  //Draw ratios
  for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
    {
      (*it).getrthshift()->SetLineColor(opts.colors[(*it).getcoli()]);
      (*it).getrthshift()->SetLineStyle(2);
      (*it).getrthshift()->SetLineWidth(2);

      (*it).getrth()->SetLineColor(opts.colors[(*it).getcoli()]);
      (*it).getrth()->SetLineWidth(2);

      vector <range> rthranges = historanges((*it).getrthshift());
      if (!opts.threepanels)
	{
	  for (vector<range>::iterator r = rthranges.begin(); r != rthranges.end(); r++)
	    {
	      (*it).getrthshift()->SetAxisRange((*r).lowedge, (*r).upedge);
	      if (!opts.onlytheory)
		(*it).Draw((TH1F*)(*it).getrthshift()->Clone(), "LX same");
	    }
	  (*it).getrthshift()->SetAxisRange((*it).getxmin(), (*it).getxmax());
	}
      
      if (!opts.points || (*it).bincenter()) //plot as continous line with dashed error bands
	{
	  for (vector<range>::iterator r = rthranges.begin(); r != rthranges.end(); r++)
	    {
	      (*it).getrth()->SetAxisRange((*r).lowedge, (*r).upedge);
	      (*it).Draw((TH1F*)(*it).getrth()->Clone(), "LX same");
	    }
	  (*it).getrth()->SetAxisRange((*it).getxmin(), (*it).getxmax());
	  if (opts.therr)
	    {
	      (*it).getrtherr()->SetLineColor(opts.colors[(*it).getcoli()]);
	      (*it).getrtherr()->SetMarkerSize(0);
	      (*it).getrtherr()->SetFillColor(opts.colors[(*it).getcoli()]);
	      (*it).getrtherr()->SetFillStyle(opts.styles[(*it).getcoli()]);
	      float toterr = 0;
	      for (int b = 1; b <= (*it).gettherr()->GetNbinsX(); b++)
		toterr += (*it).gettherr()->GetBinError(b);
	      if (toterr > 0)
		{
		  for (vector<range>::iterator r = rthranges.begin(); r != rthranges.end(); r++)
		    {
		      (*it).getrtherr()->SetAxisRange((*r).lowedge, (*r).upedge);
		      (*it).Draw((TH1F*)(*it).getrtherr()->Clone(), "E3L same");
		    }
		  (*it).getrtherr()->SetAxisRange((*it).getxmin(), (*it).getxmax());
		}
	    }
	}
      else //plot as displaced TGraphs
	{
	  TGraphAsymmErrors * r_gtherr = new TGraphAsymmErrors((*it).getrth());
	  r_gtherr->SetMarkerStyle(opts.markers[(*it).getcoli()]);
	  r_gtherr->SetLineColor(opts.colors[(*it).getcoli()]);
	  r_gtherr->SetMarkerSize(2 * opts.resolution / 1200);
	  r_gtherr->SetMarkerColor(opts.colors[(*it).getcoli()]);
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
	  r_gtherr->Draw("P same");
	}
    }	  
  
  //draw theory error borders
  if (opts.therr)
    for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
      {
	(*it).getrtherrup()->SetLineColor(opts.colors[(*it).getcoli()]);
	(*it).getrtherrdown()->SetLineColor(opts.colors[(*it).getcoli()]);
	(*it).getrtherrup()->SetLineWidth(2);
	(*it).getrtherrdown()->SetLineWidth(2);
	if (!opts.points || (*it).bincenter())
	  {
	    vector <range> rthranges = historanges((*it).getth());
	    for (vector<range>::iterator r = rthranges.begin(); r != rthranges.end(); r++)
	      {
		(*it).getrtherrup()->SetAxisRange((*r).lowedge, (*r).upedge);
		(*it).Draw((TH1F*)(*it).getrtherrup()->Clone(), "LX same");
		(*it).getrtherrdown()->SetAxisRange((*r).lowedge, (*r).upedge);
		(*it).Draw((TH1F*)(*it).getrtherrdown()->Clone(), "LX same");
	      }
	    (*it).getrtherrup()->SetAxisRange((*it).getxmin(), (*it).getxmax());
	    (*it).getrtherrdown()->SetAxisRange((*it).getxmin(), (*it).getxmax());
	  }
      }

  //Theory+shifts/Data ratio Pad (optional)
  if (opts.threepanels)
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
      for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
	mx = max(mx, (float)((*it).getrthshift()->GetMaximum()));

      mn = mx;
      for (int b = 1; b <= r_dataerr->GetNbinsX(); b++)
	r_dataerr->SetBinContent(b, r_data->GetBinContent(b) - r_data->GetBinError(b));
      if (!opts.onlytheory)
	mn = hmin(r_dataerr);
      for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
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

      //plot data
      r_templ->DrawCopy("AXIS");

      if (!opts.onlytheory)
	datahistos[0].Draw(r_data, "PE1 same");

      //plot lines at 1
      TLine *rs_one = new TLine(r_templ->GetBinLowEdge(r_templ->GetXaxis()->GetFirst()), 1, r_templ->GetXaxis()->GetBinUpEdge(r_templ->GetXaxis()->GetLast()), 1);
      rs_one->SetLineStyle(2);
      rs_one->SetLineStyle(1);
      rs_one->Draw();

      //Draw ratios
      for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
	{
	  vector <range> rthranges = historanges((*it).getrthshift());
	  for (vector<range>::iterator r = rthranges.begin(); r != rthranges.end(); r++)
	    {
	      (*it).getrthshift()->SetAxisRange((*r).lowedge, (*r).upedge);
	      (*it).Draw((TH1F*)(*it).getrthshift()->Clone(), "LX same");
	    }
	  (*it).getrthshift()->SetAxisRange((*it).getxmin(), (*it).getxmax());
	}	  
    }

  //Theory-Data pulls pad
  if (opts.twopanels || opts.threepanels)
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
      for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
	{
	  if (datahistos.size() == 1)
	    {
	      (*it).getpull()->SetFillColor(opts.colors[(*it).getcoli()]);
	      (*it).getpull()->SetFillStyle(1001);
	    }
	  (*it).getpull()->SetLineStyle(1);
	  (*it).getpull()->SetLineWidth(2);
	  (*it).getpull()->SetLineColor(opts.colors[(*it).getcoli()]);
	  datahistos[0].Draw((TH1F*)(*it).getpull()->Clone(), "same ][");
	}	  
      /*
      //redraw lines over fill area
      for (vector <dataseth>::iterator it = datahistos.begin(); it != datahistos.end(); it++)
	{
	  TH1F * redrawpull = (TH1F*)(*it).getpull()->Clone();
	  redrawpull->SetFillStyle(0);
	  redrawpull->Draw("same ][");
	}
      */
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
