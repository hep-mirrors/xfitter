#include "Chi2scanPainter.h"
#include "Outdir.h"
#include "CommandParser.h"
#include "DrawLogo.h"
#include "pdferrors.h"

#include <TF1.h>
#include <TLegend.h>

#include <vector>

//canvas of central chi2scan
vector<TCanvas*> Chi2scanPainter()
{
  vector <TCanvas*> cnvs;

  vector <string> list = chi2scanlist();
  if (list.size() == 0)
    return cnvs;

  //Plot of Chi2 scan on central PDF
  for (vector<string>::iterator it = list.begin(); it != list.end(); it++)//loop on chi2 scan labels
    {
      //Set Canvas name and margins
      char cnvname[100];
      sprintf(cnvname, "chi2scan_%s",  (*it).c_str());
      TCanvas* cnv = new TCanvas(cnvname, "", 0, 0, opts.resolution, opts.resolution);
      cnv->SetLeftMargin(lmarg);
      cnv->SetRightMargin(rmarg);
      cnv->SetTopMargin(tmarg);
      cnv->SetBottomMargin(bmarg);

      //Initialise legends
      TLegend* leg1 = new TLegend(lmarg+0.01, bmarg+0.015, lmarg+0.18,  bmarg+0.015);
      leg1->SetFillColor(0);
      leg1->SetBorderSize(0);
      leg1->SetTextAlign(12);
      leg1->SetTextSize(txtsize * 0.8);
      leg1->SetTextFont(62);
      TLegend* leg2 = new TLegend(lmarg+0.48,  bmarg+0.015, lmarg+0.48+0.15, bmarg+0.015);
      leg2->SetFillColor(0);
      leg2->SetBorderSize(0);
      leg2->SetTextAlign(12);
      leg2->SetTextSize(txtsize * 0.8);
      leg2->SetTextFont(62);

      //loop on directories
      vector <TGraph*> chi2g;
      vector <TF1*> chi2f;
      for (vector<string>::iterator itl = opts.labels.begin(); itl != opts.labels.end(); itl++)
	if (chi2scanmap[*itl].label == *it)
	  {
	    //Symmetric (pol2) or asymmetric (pol4) error treatment
	    if (!outdirs[*itl].IsAsym())
	      {
		chi2scanmap[*itl].min = chi2scanmap[*itl].min2;
		chi2scanmap[*itl].delta = chi2scanmap[*itl].delta2;
		chi2scanmap[*itl].deltap = chi2scanmap[*itl].delta2;
		chi2scanmap[*itl].deltam = chi2scanmap[*itl].delta2;
		chi2scanmap[*itl].chi2min = chi2scanmap[*itl].chi2min2;
	      }
	    else
	      {
		chi2scanmap[*itl].min = chi2scanmap[*itl].min4;
		chi2scanmap[*itl].delta = 0;
		chi2scanmap[*itl].deltap = chi2scanmap[*itl].deltap4;
		chi2scanmap[*itl].deltam = chi2scanmap[*itl].deltam4;
		chi2scanmap[*itl].chi2min = chi2scanmap[*itl].chi2min4;
	      }
      	    
	    TGraph *chi2graph = new TGraph(chi2scanmap[*itl].chi2.size());
	    int i = 0;
	    for (map<double, double>::iterator it = chi2scanmap[*itl].chi2.begin(); it != chi2scanmap[*itl].chi2.end(); it++, i++)
	      chi2graph->SetPoint(i, it->first, it->second);
	    chi2g.push_back(chi2graph);
	    chi2graph->SetMarkerStyle(opts.markers[*itl]);
	    chi2graph->SetMarkerSize(2 * opts.resolution / 1200);
	    chi2graph->SetMarkerColor(opts.colors[*itl]);

	    TF1 *parfit;
	    if (outdirs[*itl].IsAsym())
	      {
		/*
		parfit = new TF1("ParFit", "pol3");
		parfit->SetParameter(0,chi2scanmap[*itl].D);
		parfit->SetParameter(1,chi2scanmap[*itl].C);
		parfit->SetParameter(2,chi2scanmap[*itl].B);
		parfit->SetParameter(3,chi2scanmap[*itl].A);
		*/
		parfit = new TF1("ParFit", "pol4");
		parfit->SetParameter(0,chi2scanmap[*itl].e4);
		parfit->SetParameter(1,chi2scanmap[*itl].d4);
		parfit->SetParameter(2,chi2scanmap[*itl].c4);
		parfit->SetParameter(3,chi2scanmap[*itl].b4);
		parfit->SetParameter(4,chi2scanmap[*itl].a4);
	      }
	    else
	      {
		parfit = new TF1("ParFit", "[0]+(x-[2])**2/[1]**2");
		parfit->SetParameter(0,chi2scanmap[*itl].chi2min2);
		parfit->SetParameter(1,chi2scanmap[*itl].delta2);
		parfit->SetParameter(2,chi2scanmap[*itl].min2);
	      }
	    chi2f.push_back(parfit);
	    parfit->SetLineColor(opts.colors[*itl]);

	    leg1->AddEntry(chi2graph, (*itl).c_str(), "p");
	    char res[100];
	    if (outdirs[*itl].IsAsym())
	      sprintf (res, "%s + %s - %s", Round(chi2scanmap[*itl].min, min(chi2scanmap[*itl].deltap,chi2scanmap[*itl].deltam))[0].c_str(),
		       Round(chi2scanmap[*itl].deltap, min(chi2scanmap[*itl].deltap,chi2scanmap[*itl].deltam))[0].c_str(),
		       Round(chi2scanmap[*itl].deltam, min(chi2scanmap[*itl].deltap,chi2scanmap[*itl].deltam))[0].c_str());
	    else
	      sprintf (res, "%s #pm %s", Round(chi2scanmap[*itl].min, chi2scanmap[*itl].delta)[0].c_str(), Round(chi2scanmap[*itl].min, chi2scanmap[*itl].delta)[1].c_str());
	    leg2->AddEntry(parfit, res, "l");
	  }
      leg1->SetY2(bmarg+0.03+0.045*chi2g.size());
      leg2->SetY2(bmarg+0.03+0.045*chi2g.size());

      //Compute maximum an minimum for y axis
      double ymx = 0;
      double xmin, xmax, ymin, ymax;
      for (vector<TGraph*>::iterator git = chi2g.begin(); git != chi2g.end(); git++)
	{
	  (*git)->ComputeRange(xmin,ymin,xmax,ymax);
	  ymx = max(ymx,ymax);
	}
      double ymn = ymx;
      for (vector<TGraph*>::iterator git = chi2g.begin(); git != chi2g.end(); git++)
	{
	  (*git)->ComputeRange(xmin,ymin,xmax,ymax);
	  ymn = min(ymin,ymn);
	}

      double delta = ymx - ymn;
      ymx = ymx + delta*0.2;
      ymn = ymn - delta*(0.05+0.1*chi2g.size());

      //Set graphic options
      vector <TF1*>::iterator fit = chi2f.begin();
      vector<TGraph*>::iterator git = chi2g.begin();
      for (; git != chi2g.end() && fit != chi2f.end(); git++,fit++)
	{
	  if (chi2f.size() == 1)
	    {
	      (*git)->SetMarkerStyle(20);
	      (*git)->SetMarkerSize(2 * opts.resolution / 1200);
	      (*git)->SetMarkerColor(kBlack);
	      (*fit)->SetLineColor(kBlue+1);
	    }
	  
	  (*git)->SetTitle("");

	  (*fit)->SetLineStyle(2);
	  (*fit)->SetLineWidth(opts.lwidth);
	  double ymin, ymax;
	  (*git)->ComputeRange(xmin,ymin,xmax,ymax);
	  (*fit)->SetRange(xmin, xmax);
	}

      //Make template for axis
      double deltax = xmax - xmin;
      xmin = xmin-deltax*0.1;
      xmax = xmax+deltax*0.1;
      TH1F *templ = new TH1F(((string) "templ_" + cnvname).c_str(), "", 100, xmin, xmax);
      templ->GetYaxis()->SetLabelFont(62);
      templ->GetYaxis()->SetTitleFont(62);
      templ->GetYaxis()->SetLabelSize(txtsize);
      templ->GetYaxis()->SetTitleSize(txtsize);
      templ->GetYaxis()->SetTitleOffset(offset+0.3);
      templ->GetYaxis()->SetTitle("#chi^{2}");
      templ->GetXaxis()->SetLabelFont(62);
      templ->GetXaxis()->SetTitleFont(62);
      templ->GetXaxis()->SetLabelSize(txtsize);
      templ->GetXaxis()->SetTitleSize(txtsize);
      templ->GetXaxis()->SetTitle((*it).c_str());
      templ->SetStats(0);
      templ->Draw("AXIS");
      templ->SetMaximum(ymx);
      templ->SetMinimum(ymn);

      //Draw
      fit = chi2f.begin();
      git = chi2g.begin();
      for (; git != chi2g.end() && fit != chi2f.end(); git++,fit++)
	{
	  (*git)->Draw("P same");
	  (*fit)->Draw("same");
	}
      leg1->Draw();
      leg2->Draw();
      DrawLabels();
      if (opts.drawlogo)
        DrawLogo()->Draw();   

      //Store Canvas
      cnvs.push_back(cnv);
    }

  //Make plots of PDF uncertainties
  for (vector<string>::iterator it = list.begin(); it != list.end(); it++)//loop on chi2 scan labels
    {
      //Set Canvas name and margins
      char cnvname[100];
      sprintf(cnvname, "deltapdf_%s",  (*it).c_str());
      TCanvas* cnv = new TCanvas(cnvname, "", 0, 0, opts.resolution, opts.resolution);
      cnv->SetLeftMargin(lmarg+0.01);
      cnv->SetRightMargin(rmarg);
      cnv->SetTopMargin(tmarg);
      cnv->SetBottomMargin(bmarg);

      //Initialise legends
      TLegend* leg1 = new TLegend(lmarg+0.018, bmarg+0.015, lmarg+0.15,  bmarg+0.015);
      leg1->SetFillColor(0);
      leg1->SetBorderSize(0);
      leg1->SetTextAlign(12);
      leg1->SetTextSize(txtsize * 0.8);
      leg1->SetTextFont(62);
      TPaveText* leg2 = new TPaveText(lmarg+0.48,  bmarg+0.015, 1-rmarg-0.01, bmarg+0.015, "NDC");
      leg2->SetFillColor(0);
      leg2->SetBorderSize(0);
      leg2->SetTextAlign(12);
      leg2->SetTextSize(txtsize * 0.8);
      leg2->SetTextFont(62);

      //loop on directories
      vector <TH1F*> pdfdelta;
      for (vector<string>::iterator itl = opts.labels.begin(); itl != opts.labels.end(); itl++)
	if (chi2scanmap[*itl].label == *it)
	  {
	    if (chi2scanmap[*itl].pdfchi2.size() == 0)
	      continue;

	    for (unsigned int i = 1; i <= chi2scanmap[*itl].pdfmin2.size(); i++)
	      {
		if (!outdirs[*itl].IsAsym())
		  chi2scanmap[*itl].pdfmin[i] = chi2scanmap[*itl].pdfmin2[i];
		else
		  chi2scanmap[*itl].pdfmin[i] = chi2scanmap[*itl].pdfmin4[i];
	      }
	    
	    TH1F* h_delta = new TH1F(((string)"h_delta_" + *itl).c_str(),"", chi2scanmap[*itl].pdfchi2.size(), 1, chi2scanmap[*itl].pdfchi2.size()+1);
	    vector <double> xi;
	    // no need to distinguish between pol2 and pol3-4 fits here //if (outdirs[*itl].IsAsym())
	    xi.push_back(chi2scanmap[*itl].min);
	    for (unsigned int i = 1; i <= chi2scanmap[*itl].pdfmin.size(); i++)
	      {
		h_delta->SetBinContent(i, chi2scanmap[*itl].pdfmin[i]);
		xi.push_back(chi2scanmap[*itl].pdfmin[i]);
		char member[10];
		sprintf (member, "%d", i);
		h_delta->GetXaxis()->SetBinLabel(i, member);
	      }

	    pdfdelta.push_back(h_delta);
	    h_delta->SetLineColor(opts.colors[*itl]);
	    if (opts.bw)
	      h_delta->SetLineStyle(opts.lstyles[*itl]);


	    leg1->AddEntry(h_delta, (*itl).c_str(), "l");

	    double eplus, eminus;
	    //Ugly patch do distinguish between MC and hessian PDFs (-> assume they are MC for N > 80)
	    if (xi.size() > 80)
	      //deltaasym(xi, median(xi), eplus, eminus, cl(1));
	      eplus = eminus = rms(xi);
	    else if (xi.size() == 42) //Uglier patch for HERAPDF2.0
	      vardeltaasym(xi, 3, eplus, eminus);
	    else
	      ahessdeltaasym(xi, eplus, eminus);

	    char res[100];
	    double mn;
	    // no need to distinguish between pol2 and pol3-4 fits here //if (outdirs[*itl].IsAsym())
	    mn = chi2scanmap[*itl].min;
	    
	    sprintf (res, "%s +%s -%s", 
		     Round(mn, min(eplus,eminus))[0].c_str(), 
		     Round(eplus, min(eplus,eminus))[0].c_str(), 
		     Round(eminus, min(eplus,eminus))[0].c_str());
	    leg2->AddText(res);
	  }
      if (pdfdelta.size() == 0)
	continue;
      leg1->SetY2(bmarg+0.03+0.045*pdfdelta.size());
      leg2->SetY2(bmarg+0.03+0.045*pdfdelta.size());

      //Compute maximum and minimum for y axis
      double ymx = (*pdfdelta.begin())->GetMaximum();
      double ymn = (*pdfdelta.begin())->GetMinimum();
      for (vector<TH1F*>::iterator it = pdfdelta.begin(); it != pdfdelta.end(); it++)
	{
	  ymx = max((*it)->GetMaximum(),ymx);
	  ymn = min((*it)->GetMinimum(),ymn);
	}
      double delta = ymx - ymn;
      ymx = ymx + delta*0.25;
      ymn = ymn - delta*(0.05+0.11*pdfdelta.size());

      //Compute maximum and minimum for x axis
      double xmx = (*pdfdelta.begin())->GetXaxis()->GetBinUpEdge((*pdfdelta.begin())->GetXaxis()->GetLast());
      double xmn = (*pdfdelta.begin())->GetXaxis()->GetBinLowEdge((*pdfdelta.begin())->GetXaxis()->GetFirst());
      for (vector<TH1F*>::iterator it = pdfdelta.begin(); it != pdfdelta.end(); it++)
	{
	  xmx = max((*it)->GetXaxis()->GetBinUpEdge((*it)->GetXaxis()->GetLast()), xmx);
	  xmn = min((*it)->GetXaxis()->GetBinLowEdge((*it)->GetXaxis()->GetFirst()), xmn);
	}

      //Set graphic options
      for (vector<TH1F*>::iterator it = pdfdelta.begin(); it != pdfdelta.end(); it++)
	{
	  if (pdfdelta.size() == 1)
	    (*it)->SetLineColor(kBlack);

	  (*it)->SetStats(0);
	  (*it)->SetLineWidth(opts.lwidth);
	  (*it)->SetTitle("");
	}

      //Make template for axis
      TH1F *templ = new TH1F(((string) "templ-pdf_" + cnvname).c_str(), "", xmx-xmn, xmn, xmx);
      for (unsigned int i = 1; i < xmx; i++)
	{
	  if (i % 5 == 1)
	    {
	      char member[10];
	      sprintf (member, "%d", i);
	      templ->GetXaxis()->SetBinLabel(i, member);
	    }
	  else
	      templ->GetXaxis()->SetBinLabel(i, "");
	}
      templ->GetXaxis()->CenterLabels();
      templ->GetXaxis()->LabelsOption("v");
      templ->GetYaxis()->SetLabelFont(62);
      templ->GetYaxis()->SetTitleFont(62);
      templ->GetYaxis()->SetLabelSize(txtsize);
      templ->GetYaxis()->SetTitleSize(txtsize);
      templ->GetYaxis()->SetTitleOffset(offset+0.4);
      templ->GetYaxis()->SetTitle((*it).c_str());
      templ->GetXaxis()->SetLabelFont(62);
      templ->GetXaxis()->SetTitleFont(62);
      templ->GetXaxis()->SetLabelSize(txtsize*1.5);
      templ->GetXaxis()->SetTitleSize(txtsize);
      //templ->GetXaxis()->SetTitle((lhapdfset + " member").c_str());
      templ->GetXaxis()->SetTitle("PDF member");
      templ->SetStats(0);
      templ->Draw("AXIS");
      templ->SetMaximum(ymx);
      templ->SetMinimum(ymn);

      //Draw
      for (vector<TH1F*>::iterator it = pdfdelta.begin(); it != pdfdelta.end(); it++)
	(*it)->Draw("hist same ][");
      leg2->Draw();
      leg1->Draw();
      DrawLabels();
      if (opts.drawlogo)
        DrawLogo()->Draw();   

      //Store Canvas
      cnvs.push_back(cnv);
    }

  return cnvs;
}
