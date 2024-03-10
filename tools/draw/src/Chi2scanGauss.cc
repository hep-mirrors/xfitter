#include "Chi2scanGauss.h"
#include "Outdir.h"
#include "CommandParser.h"
#include "DrawLogo.h"
#include "pdferrors.h"

#include <TF1.h>
#include <TLegend.h>

#include <vector>

//canvas of central chi2scan
vector<TCanvas*> Chi2scanGauss()
{
  vector <TCanvas*> cnvs;

  vector <string> list = chi2scanlist();
  if (list.size() == 0)
    return cnvs;

  //Plot of Chi2 scan gauss fit of MC toys
  for (vector<string>::iterator it = list.begin(); it != list.end(); it++) //loop on chi2 scan labels (usually one)
    {
      //Set Canvas name and margins
      char cnvname[100];
      sprintf(cnvname, "chi2scan_gauss_%s",  (*it).c_str());
      TCanvas* cnv = new TCanvas(cnvname, "", 0, 0, opts.resolution, opts.resolution);
      cnv->SetLeftMargin(lmarg);
      cnv->SetRightMargin(rmarg);
      cnv->SetTopMargin(tmarg);
      cnv->SetBottomMargin(bmarg);

      
      //Initialise legends
      TLegend *leg = new TLegend(0.11, 0.75, 0.25, 0.89);
      leg->SetBorderSize(0); 
      leg->SetTextSize(0.035); 
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

      double xmin, xmax;
      
      //loop on directories
      vector <TH1D*> hgauss;
      vector <TF1*> fgauss;
      for (vector<string>::iterator itl = opts.labels.begin(); itl != opts.labels.end(); itl++)
	if (chi2scanmap[*itl].label == *it)
	  if (chi2scanmap[*itl].min2_mc.size() > 0)
	    {
	      //Symmetric (pol2) or asymmetric (pol4) error treatment
	      if (!outdirs[*itl].IsAsym())
		{
		  chi2scanmap[*itl].min = mean(chi2scanmap[*itl].min2_mc);
		  chi2scanmap[*itl].delta = mean(chi2scanmap[*itl].delta2_mc);
		  chi2scanmap[*itl].deltap = mean(chi2scanmap[*itl].delta2_mc);
		  chi2scanmap[*itl].deltam = mean(chi2scanmap[*itl].delta2_mc);
		  //		  chi2scanmap[*itl].chi2min = mean(chi2scanmap[*itl].chi2min2_mc);
		}
	      else
	      {
		chi2scanmap[*itl].min = mean(chi2scanmap[*itl].min4_mc);
		chi2scanmap[*itl].delta = 0;
		chi2scanmap[*itl].deltap = mean(chi2scanmap[*itl].deltap4_mc);
		chi2scanmap[*itl].deltam = mean(chi2scanmap[*itl].deltam4_mc);
		//		chi2scanmap[*itl].chi2min = mean(chi2scanmap[*itl].chi2min4_mc);
	      }

	      
	      map <double, double>::iterator cf = chi2scanmap[*itl].chi2.begin();
	      map <double, double>::iterator cl = chi2scanmap[*itl].chi2.end(); cl--;
	      xmin = cf->first;
	      xmax = cl->first;
	      TH1D * hg = new TH1D(((string) "gauss_" + cnvname + (*itl)).c_str(), "", 3*sqrt(chi2scanmap[*itl].min2_mc.size()), xmin, xmax);
	      //loop on mc replicas
	      for (int i = 0; i < chi2scanmap[*itl].min2_mc.size(); i++)
		hg->Fill(chi2scanmap[*itl].min4_mc[i]);


	      TF1 *fg = new TF1(((string) "gauss_func_" + cnvname + (*itl)).c_str(), "gaus", xmin, xmax);
	      hg->Fit(((string) "gauss_func_" + cnvname + (*itl)).c_str(), "WLN");
	      
	      hg->SetLineColor(opts.colors[*itl]);
	      fg->SetLineColor(opts.colors[*itl]);
	      hg->SetXTitle(chi2scanmap[*itl].label.c_str());
	      hg->SetYTitle("MC toys");
	      hg->SetStats(0);

	      hgauss.push_back(hg);
	      fgauss.push_back(fg);


	      char text[200];
	      //sprintf(text, "%s %.2f #pm %.2f", (*itl).c_str(), fg->GetParameter(1), fg->GetParameter(2));
	      sprintf(text, "%s %s #pm %s", Round(fg->GetParameter(1), fg->GetParameter(2))[0].c_str(), Round(fg->GetParameter(1), fg->GetParameter(2))[1].c_str()," ");
	      leg->AddEntry(hg, text);

	      leg1->AddEntry(hg, (*itl).c_str(), "p");
	      leg2->AddEntry(fg, text, "l");
	      
	    } //end loop on directories

      if (hgauss.size() == 0) //No MC replica sets of dir found
	continue;
      
      //Compute maximum an minimum for y axis
      double ymax = 0;
      double ymin = 0;
      for (vector<TH1D*>::iterator hit = hgauss.begin(); hit != hgauss.end(); hit++)
	{
	  double ymx = (*hit)->GetMaximum();
	  ymax = max(ymx,ymax);
	}
      /*
      double ymin = ymax;
      for (vector<TGraph*>::iterator hit = hgauss.begin(); hit != hgauss.end(); hit++)
	{
	  double ymn = (*hit)->GetMinimum();
	  ymin = min(ymin,ymn);
	}
      */
      double delta = ymax - ymin;
      ymax = ymax + delta*0.2;
      
      /*
      //lines...
 TLine *line1 = new TLine(f1->GetParameter(1),0,f1->GetParameter(1),ymax);
 line1->SetLineColor(kRed);
 line1->SetLineWidth(2);
 line1->SetLineStyle(2);
 line1->Draw("same");
      */
	      
      //Set graphic options
      vector <TF1*>::iterator fit = fgauss.begin();
      vector <TH1D*>::iterator hit = hgauss.begin();
      for (; hit != hgauss.end() && fit != fgauss.end(); hit++,fit++)
	{
	  (*hit)->SetTitle("");

	  (*fit)->SetLineStyle(2);
	  (*fit)->SetLineWidth(opts.lwidth);
	  //(*hit)->SetMaximum(ymax);
	  //(*hit)->SetMinimum(ymin);
	}

      //leg1->SetY2(bmarg+0.03+0.045*chi2g.size());
      //leg2->SetY2(bmarg+0.03+0.045*chi2g.size());

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
      templ->GetYaxis()->SetTitle("MC toys");
      templ->GetXaxis()->SetLabelFont(62);
      templ->GetXaxis()->SetTitleFont(62);
      templ->GetXaxis()->SetLabelSize(txtsize);
      templ->GetXaxis()->SetTitleSize(txtsize);
      templ->GetXaxis()->SetTitle((*it).c_str());
      templ->SetStats(0);
      templ->Draw("AXIS");
      templ->SetMaximum(ymax);
      templ->SetMinimum(ymin);

      //Draw
      fit = fgauss.begin();
      hit = hgauss.begin();
      for (; hit != hgauss.end() && fit != fgauss.end(); hit++,fit++)
	{
	  (*hit)->Draw("e0 same");
	  (*fit)->Draw("same");
	}
 
      leg->Draw();
      //leg1->Draw();
      //leg2->Draw();
      DrawLabels();
      if (opts.drawlogo)
        DrawLogo()->Draw();   

      //Store Canvas
      cnvs.push_back(cnv);
    }

  return cnvs;
}
