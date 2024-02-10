#include "PdfsPainter.h"
#include "DrawLogo.h"
#include "CommandParser.h"
#include "Outdir.h"

#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TMath.h>

#include <iostream>


vector <TCanvas*> PdfsPainter(double q2, pdftype ipdf)
{
  vector <TCanvas*> cnvs;

  char q2str[30];
  if (q2 < 10)
    sprintf(q2str, "%.1f",  q2);
  else
    sprintf(q2str, "%.0f",  q2);

  vector <TGraphAsymmErrors*> pdfgraphs;
  vector <string> labels;
  TGraphAsymmErrors* pdfempunc[6][3];
  bool blIs3bands[6];
  for (vector<string>::iterator itl = opts.labels.begin(); itl != opts.labels.end(); itl++)
    if (pdfmap[*itl].Central.find(q2) !=  pdfmap[*itl].Central.end())
      {
        char pdfname[80];
        blIs3bands[itl - opts.labels.begin()] = outdirs[*itl].Is3bands();
        if (!blIs3bands[itl - opts.labels.begin()]) {
          pdfgraphs.push_back(pdfmap[*itl].Central[q2].GetPdf(ipdf));
          labels.push_back(*itl);

          sprintf(pdfname, "dir%ld_q2_%s_pdf_%s", (itl-opts.labels.begin()+1), q2str, pdffiles[ipdf].c_str());
          pdfgraphs.back()->SetName(pdfname);
        } else {
          pdfgraphs.push_back(pdfmap[*itl].Central[q2].GetPdfCen(ipdf));
          labels.push_back(*itl);
          sprintf(pdfname, "dir%ld_q2_%s_pdf_%s_cen", (itl-opts.labels.begin()+1), q2str, pdffiles[ipdf].c_str());
          pdfgraphs.back()->SetName(pdfname);

          pdfempunc[itl - opts.labels.begin()][0] = pdfmap[*itl].Central[q2].GetPdfParam(ipdf);
          sprintf(pdfname, "dir%ld_q2_%s_pdf_%s_param", (itl-opts.labels.begin()+1), q2str, pdffiles[ipdf].c_str());
          pdfempunc[itl - opts.labels.begin()][0]->SetName(pdfname);

          pdfempunc[itl - opts.labels.begin()][1] = pdfmap[*itl].Central[q2].GetPdfModel(ipdf);
          sprintf(pdfname, "dir%ld_q2_%s_pdf_%s_model", (itl-opts.labels.begin()+1), q2str, pdffiles[ipdf].c_str());
          pdfempunc[itl - opts.labels.begin()][1]->SetName(pdfname);

          pdfempunc[itl - opts.labels.begin()][2] = pdfmap[*itl].Central[q2].GetPdfExp(ipdf);
          sprintf(pdfname, "dir%ld_q2_%s_pdf_%s_exp", (itl-opts.labels.begin()+1), q2str, pdffiles[ipdf].c_str());
          pdfempunc[itl - opts.labels.begin()][2]->SetName(pdfname);
        }
      }

  for (vector <TGraphAsymmErrors*>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
      allgraphs.push_back(*it);

  if (pdfgraphs.size() < 1)
    {
      cout << "Error: Empty pdf TGraph vector for pdf: " << pdffiles[(int)ipdf] << endl;
      exit(1);
    }

  char cnvname[60];
  sprintf(cnvname, "q2_%s_pdf_%s",  q2str, pdffiles[ipdf].c_str());

  if(opts.xmin==-1&&opts.xmax==-1){
    //Set xmin xmax
    opts.xmin= 1e100;
    opts.xmax=-1e100;
    for(const auto&graph:pdfgraphs){
      if(!graph)continue;
      TAxis*ax=graph->GetXaxis();
      double xmin = ax->GetXmin();
      if(xmin < 0.0 || (xmin == 0.0 && opts.logx))
        xmin = graph->GetX()[0];
      opts.xmin = min(opts.xmin, xmin);
      double xmax = ax->GetXmax();
      if(xmax > 1.0)
        xmax = graph->GetX()[graph->GetN() - 1];
      opts.xmax=max(opts.xmax, xmax);
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
              if (xi < opts.xmin || xi > opts.xmax)
                {
                  (*it)->RemovePoint(i);
                  if (blIs3bands[it - pdfgraphs.begin()]) {
                    pdfempunc[it - pdfgraphs.begin()][0]->RemovePoint(i);
                    pdfempunc[it - pdfgraphs.begin()][1]->RemovePoint(i);
                    pdfempunc[it - pdfgraphs.begin()][2]->RemovePoint(i);
                  }
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

      if (blIs3bands[it - pdfgraphs.begin()]) {
        pdfempunc[it - pdfgraphs.begin()][0]->SetFillColor(3);
        pdfempunc[it - pdfgraphs.begin()][1]->SetFillColor(5);
        pdfempunc[it - pdfgraphs.begin()][2]->SetFillColor(2);
        for (int iun=0; iun<3; iun++) pdfempunc[it - pdfgraphs.begin()][iun]->SetFillStyle(1001);
      }

      (*it)->SetFillColor(opts.colors[labels[it-pdfgraphs.begin()]]);
      if (opts.filledbands || opts. transparentbands)
	{
	  (*it)->SetFillStyle(1001);
	  if (opts. transparentbands)
	    (*it)->SetFillColorAlpha((*it)->GetFillColor(), 0.5);
	}
      else
        (*it)->SetFillStyle(opts.styles[labels[it-pdfgraphs.begin()]]);
      (*it)->SetLineStyle(1);
      (*it)->SetLineWidth(opts.lwidth);
      (*it)->SetLineColor(opts.colors[labels[it-pdfgraphs.begin()]]);
      colindx++;
    }

  //Calculate maximum and minimum of y axis
  double mx = 0;
  for (vector <TGraphAsymmErrors*>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    for (int i = 0; i < (*it)->GetN(); i++)
      {
        double xi = (*it)->GetX()[i];
        double val = (*it)->GetY()[i];
        double errhigh;
        if (!blIs3bands[it - pdfgraphs.begin()]) {
          errhigh = (*it)->GetErrorYhigh(i);
        } else {
          errhigh = pdfempunc[it - pdfgraphs.begin()][0]->GetErrorYhigh(i);
        }
        if (xi >= opts.xmin && xi <= opts.xmax)
          mx = max(mx, val+errhigh);
      }
  double mn = mx;
  for (vector <TGraphAsymmErrors*>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    for (int i = 0; i < (*it)->GetN(); i++)
      {
        double xi = (*it)->GetX()[i];
        double val = (*it)->GetY()[i];
        double errlow;
        if (!blIs3bands[it - pdfgraphs.begin()]) {
          errlow = (*it)->GetErrorYlow(i);
        } else {
          errlow = pdfempunc[it - pdfgraphs.begin()][0]->GetErrorYlow(i);
        }
        if (xi >= opts.xmin && xi <= opts.xmax)
          mn = min(mn, val-errlow);
      }

  //Prepare TGraphs
  TMultiGraph * mg = new TMultiGraph(((string)cnvname + "_multigraph").c_str(), "");
  TMultiGraph * mg_lines = new TMultiGraph(((string)cnvname + "_multigraph_lines").c_str(), "");
  TMultiGraph * mg_dotted_lines = new TMultiGraph(((string)cnvname + "_multigraph_dotted_lines").c_str(), "");
  TMultiGraph * mg_shade = new TMultiGraph(((string)cnvname + "_multigraph_shade").c_str(), "");
  for (vector <TGraphAsymmErrors*>::iterator it = pdfgraphs.begin(); it != pdfgraphs.end(); it++)
    {
      if (*it == 0)
        continue;
      //Prepare graph line borders and graph shade
      int npoints = (*it)->GetN();
      double val_x[npoints], val_y[npoints], val_high_y[npoints], val_low_y[npoints];
      double xsh[2*npoints], ysh[2*npoints], yshEMP[3][2*npoints];

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

          if (blIs3bands[it - pdfgraphs.begin()]) {
            for (int iun=0; iun<3; iun++) {
              yshEMP[iun][i] = val + pdfempunc[it - pdfgraphs.begin()][iun]->GetErrorYhigh(i);
              yshEMP[iun][npoints + i] = (*it)->GetY()[npoints-i-1] - pdfempunc[it - pdfgraphs.begin()][iun]->GetErrorYlow(npoints-i-1);
            }
          }
        }

      TGraph *centr = new TGraph(npoints, val_x, val_y);
      TGraph *high = new TGraph(npoints, val_x, val_high_y);
      TGraph *low = new TGraph(npoints, val_x, val_low_y);
      TGraph *shade = new TGraph(2*npoints, xsh, ysh);

      TGraph *high_dot = new TGraph(npoints, val_x, val_high_y);
      TGraph *low_dot = new TGraph(npoints, val_x, val_low_y);

      TGraph *high_shade = new TGraph(npoints, val_x, val_high_y);
      TGraph *low_shade = new TGraph(npoints, val_x, val_low_y);

      TGraph *shadeEMP[3];
      if (blIs3bands[it - pdfgraphs.begin()]) {
        for (int iun=0; iun<3; iun++) {
          shadeEMP[iun] = new TGraph(2*npoints, xsh, yshEMP[iun]);
          shadeEMP[iun]->SetFillColor(pdfempunc[it - pdfgraphs.begin()][iun]->GetFillColor());
          shadeEMP[iun]->SetFillStyle(pdfempunc[it - pdfgraphs.begin()][iun]->GetFillStyle());
          shadeEMP[iun]->SetLineWidth(0);
        }
      }

      //Set border lines and shade fill
      centr->SetLineColor((*it)->GetLineColor());
      centr->SetLineStyle(1);
      centr->SetLineWidth(opts.lwidth);
      high->SetLineColor((*it)->GetLineColor());
      high->SetLineWidth(opts.lwidth);
      high->SetLineStyle(1);
      high_dot->SetLineWidth(1);
      high_dot->SetLineColor((*it)->GetLineColor());
      high_dot->SetLineStyle(2);
      high_shade->SetLineWidth(1);
      high_shade->SetLineColor((*it)->GetLineColor());
      high_shade->SetLineStyle(1);
      low->SetLineColor((*it)->GetLineColor());
      low->SetLineWidth(opts.lwidth);
      low->SetLineStyle(1);
      low_dot->SetLineWidth(1);
      low_dot->SetLineColor((*it)->GetLineColor());
      low_dot->SetLineStyle(2);
      low_shade->SetLineWidth(1);
      low_shade->SetLineColor((*it)->GetLineColor());
      low_shade->SetLineStyle(1);
      shade->SetLineColor((*it)->GetLineColor());
      shade->SetFillColor((*it)->GetLineColor());
      shade->SetFillStyle((*it)->GetFillStyle());
      if (opts. transparentbands)
	shade->SetFillColorAlpha((*it)->GetFillColor(), 0.5);
      shade->SetLineWidth(0);

      //add graphs
      mg->Add((*it));
      mg_lines->Add(centr);
      mg_lines->Add(high);
      mg_lines->Add(low);
      if (it+1 != pdfgraphs.end())
        {
          mg_dotted_lines->Add(high_dot);
          mg_dotted_lines->Add(low_dot);
        }
      mg_shade->Add(shade, "f");
      if (blIs3bands[it - pdfgraphs.begin()]) {
        for (int iun=0; iun<3; iun++) {
          mg_shade->Add(shadeEMP[iun], "f");
        }
      }
      mg_shade->Add(high_shade, "l");
      mg_shade->Add(low_shade, "l");
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
  mg->SetTitle(((string)" ; x  ; x" + pdflabels[ipdf] + "(x," + opts.q2label + ")").c_str());

  mg->Draw("A"); //need to draw with A option to create axis

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
  mg->GetXaxis()->SetTitleFont(opts.rootfont);
  mg->GetXaxis()->SetLabelFont(opts.rootfont);
  mg->GetXaxis()->SetTitleSize(txtsize);
  mg->GetXaxis()->SetLabelSize(txtsize);

  mg->GetYaxis()->SetTitleFont(opts.rootfont);
  mg->GetYaxis()->SetLabelFont(opts.rootfont);
  mg->GetYaxis()->SetTitleSize(txtsize);
  mg->GetYaxis()->SetLabelSize(txtsize);
  mg->GetYaxis()->SetTitleOffset(offset);

  mg_shade->Draw("");
  if (opts.filledbands && !opts.transparentbands)
    mg_dotted_lines->Draw("l");
  else
    mg_lines->Draw("l");

  //Make legend
  //TLegend * leg = new TLegend(lmarg+0.03, 1-tmarg-0.05-pdfgraphs.size()*0.05, lmarg+0.33, 1-tmarg-0.01);
  TLegend * leg = new TLegend(lmarg+0.03, 1-tmarg-0.07-pdfgraphs.size()*0.05, lmarg+0.33, 1-tmarg-0.03);
  leg->SetTextFont(opts.rootfont);
  leg->SetTextSize(txtsize);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  string q2string = opts.q2label + " = " + q2str + " GeV^{2}";
  if (fabs(sqrt(q2) - 80.385) < 1)
    q2string = opts.q2label + " = m_{W}^{2}";
  else if (fabs(sqrt(q2) - 91.1876) < 1)
    q2string = opts.q2label + " = m_{Z}^{2}";
  else if (fabs(sqrt(q2) - 125) < 1)
    q2string = opts.q2label + " = m_{H}^{2}";
  else if (fabs(sqrt(q2) - 173.3) < 1)
    q2string = opts.q2label + " = m_{t}^{2}";
  //  leg->AddEntry((TObject*)0, ((string)"x" + pdflabels[ipdf] + " - " + opts.q2label + " = " + q2string).c_str(), "");
  leg->AddEntry((TObject*)0, q2string.c_str(), "");

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
  if (opts.drawlogo)
    DrawLogo()->Draw();
  DrawLabels("ur");

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

  //Calculate maximum and minimum of x axis
  double xmnforbds, xmxforbds;
  for (vector <TGraphAsymmErrors*>::iterator it = rlist.begin(); it != rlist.end(); it++)
    {
      if (*it == 0)
        continue;

      TGraphAsymmErrors* r = *it;
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
    }

  //Calculate maximum and minimum of y axis
  mx = 1;
  if (opts.abserror)
    mx = 0;
  for (vector <TGraphAsymmErrors*>::iterator it = rlist.begin(); it != rlist.end(); it++)
    {
      if (*it == 0)
        continue;

      TGraphAsymmErrors* r = *it;
      if (opts.rmax == 0 && opts.rmin == 0)
        for (int i = 0; i < r->GetN(); i++)
          {
            double xi = r->GetX()[i];
            double yi_h = r->GetY()[i] + fabs(r->GetErrorYhigh(i));
            if (xi >= xmnforbds && xi <= xmxforbds)
              mx = max(mx, yi_h);
          }
      else
        mx = opts.rmax;
    }
  mn = mx;
  for (vector <TGraphAsymmErrors*>::iterator it = rlist.begin(); it != rlist.end(); it++)
    {
      if (*it == 0)
        continue;

      TGraphAsymmErrors* r = *it;
      if (opts.rmax == 0 && opts.rmin == 0)
        for (int i = 0; i < r->GetN(); i++)
          {
            double xi = r->GetX()[i];
            double yi_l = r->GetY()[i] - fabs(r->GetErrorYlow(i));
            if (xi >= xmnforbds && xi <= xmxforbds)
              mn = min(mn, yi_l);
          }
      else
        mn = opts.rmin;
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
  TMultiGraph * mg_ratio_dotted_lines = new TMultiGraph(((string)cnvname + "_multigraph_ratio_dotted_lines").c_str(), "");
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
          val_y[i] = ratio;
          val_high_y[i] = high;
          val_low_y[i] = low;
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

      TGraph *r_high_dot = new TGraph(npoints, val_x, val_high_y);
      TGraph *r_low_dot = new TGraph(npoints, val_x, val_low_y);

      TGraph *r_high_shade = new TGraph(npoints, val_x, val_high_y);
      TGraph *r_low_shade = new TGraph(npoints, val_x, val_low_y);


      //Set border lines and shade fill
      r_centr->SetLineColor(r->GetLineColor());
      r_centr->SetLineStyle(1);
      r_centr->SetLineWidth(opts.lwidth);
      r_high->SetLineColor(r->GetLineColor());
      r_high->SetLineStyle(1);
      r_high->SetLineWidth(opts.lwidth);
      r_high_dot->SetLineWidth(1);
      r_high_dot->SetLineColor(r->GetLineColor());
      r_high_dot->SetLineStyle(2);
      r_high_shade->SetLineWidth(1);
      r_high_shade->SetLineColor(r->GetLineColor());
      r_high_shade->SetLineStyle(1);
      r_low->SetLineColor(r->GetLineColor());
      r_low->SetLineStyle(1);
      r_low->SetLineWidth(opts.lwidth);
      r_low_dot->SetLineWidth(1);
      r_low_dot->SetLineColor(r->GetLineColor());
      r_low_dot->SetLineStyle(2);
      r_low_shade->SetLineWidth(1);
      r_low_shade->SetLineColor(r->GetLineColor());
      r_low_shade->SetLineStyle(1);
      r_shade->SetLineColor(r->GetLineColor());
      r_shade->SetFillColor(r->GetLineColor());
      r_shade->SetFillStyle(r->GetFillStyle());
      if (opts. transparentbands)
	r_shade->SetFillColorAlpha(r->GetFillColor(), 0.5);
      r_shade->SetLineWidth(0);


      //add graphs
      mg_ratio->Add(r);
      mg_ratio_lines->Add(r_centr);
      mg_ratio_lines->Add(r_high);
      mg_ratio_lines->Add(r_low);
      if (it+1 != pdfgraphs.end())
        {
          mg_ratio_dotted_lines->Add(r_high_dot);
          mg_ratio_dotted_lines->Add(r_low_dot);
        }
      mg_ratio_shade->Add(r_shade, "f");
      mg_ratio_shade->Add(r_high_shade, "l");
      mg_ratio_shade->Add(r_low_shade, "l");
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
  mg_ratio->Draw("A"); //Create axis
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
  mg_ratio->GetYaxis()->SetTitle(((string)" x" + pdflabels[ipdf] + "(x," + opts.q2label + ")/x" + pdflabels[ipdf] + "(x," + opts.q2label + ")_{ref}").c_str());
  if (opts.relerror)
    mg_ratio->GetYaxis()->SetTitle(((string)" #deltax" + pdflabels[ipdf] + "/x" + pdflabels[ipdf]).c_str());
  if (opts.abserror)
    mg_ratio->GetYaxis()->SetTitle(((string)" #deltax" + pdflabels[ipdf] + "").c_str());


  mg_ratio->GetXaxis()->Set(100, opts.xmin, opts.xmax);
  mg_ratio->GetXaxis()->SetTitleFont(opts.rootfont);
  mg_ratio->GetXaxis()->SetLabelFont(opts.rootfont);
  mg_ratio->GetXaxis()->SetTitleSize(txtsize);
  mg_ratio->GetXaxis()->SetLabelSize(txtsize);
  //  mg_ratio->GetXaxis()->SetTitleOffset(offset);
  
  mg_ratio->GetYaxis()->SetTitleFont(opts.rootfont);
  mg_ratio->GetYaxis()->SetLabelFont(opts.rootfont);
  mg_ratio->GetYaxis()->SetTitleSize(txtsize);
  mg_ratio->GetYaxis()->SetLabelSize(txtsize);
  mg_ratio->GetYaxis()->SetTitleOffset(offset);
  mg_ratio->GetYaxis()->SetNdivisions(506);

  //  mg_ratio->Draw("ALE3");
  mg_ratio_shade->Draw("");
  if (opts.filledbands && !opts.transparentbands)
    mg_ratio_dotted_lines->Draw("l");
  else
    mg_ratio_lines->Draw("l");

  //Make legend
  TLegend * leg2;
  if (ipdf == ubar || ipdf == dbar || ipdf == s || ipdf == Sea || ipdf == g)
    //leg2 = new TLegend(lmarg+0.03, 1-tmarg-0.05-pdfgraphs.size()*0.05, lmarg+0.33, 1-tmarg-0.01);
    leg2 = new TLegend(lmarg+0.03, 1-tmarg-0.07-pdfgraphs.size()*0.05, lmarg+0.33, 1-tmarg-0.03);
  else
    //leg2 = new TLegend(lmarg+0.18, 1-tmarg-0.05-pdfgraphs.size()*0.05, lmarg+0.45, 1-tmarg-0.01);
    leg2 = new TLegend(lmarg+0.18, 1-tmarg-0.07-pdfgraphs.size()*0.05, lmarg+0.45, 1-tmarg-0.03);
  leg2->SetTextFont(opts.rootfont);
  leg2->SetTextSize(txtsize);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->AddEntry((TObject*)0, q2string.c_str(), "");

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

  if (opts.drawlogo)
    DrawLogo("bc")->Draw();
  //DrawLabels("bc");
  DrawLabels("ur");

  return cnvs;
}
