#include "ShiftPainter.h"
#include "CommandParser.h"
#include "DrawLogo.h"
#include "Outdir.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>

#include <stdlib.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TLine.h>
#include <math.h>

#include "FileOpener.h"

struct shtype
{
  double val;
  double err;
  int dataid;
};
typedef map<string,vector<shtype> > shlisttype;

using namespace std;
vector <TCanvas*> ShiftPainter(vector<string> dirs)
{
  vector <TCanvas*> shiftscnv;

  //optimise settings
  if (opts.adjshift)
    {
      //heigth reserved for each shift in points
      opts.shgth = 30 + dirs.size() * 10;
      //number of shifts per plot
      opts.spp = (int) (1200. * (1 - tmarg - bmarg) / opts.shgth);
    }

  //make shifts list
  shlisttype shlist;
  //read shifts
  for (vector<string>::iterator itl = opts.labels.begin(); itl != opts.labels.end(); itl++)
    {
      if (outdirs[*itl].IsMCreplica())
        continue;

      //      if (itl->find("profiled")) {
      //        continue;
      //      }

      // string fname = outdirs[*itl].GetName() + "/Results.txt";
      // ifstream f(fname.c_str());
      // if (!f.good())
        // {
          // cout << "File " << fname << " is empty (or io error) (in ShiftPainter)" << endl;
          // continue;
        // }
        
      string dirname = outdirs[*itl].GetName();
      InFileOpener_t fo;
      fo.Add(dirname + "/Results.txt");
      fo.Add(dirname + "/Results_0.txt");
      if(fo.Open()) continue;
      ifstream &f = fo.GetStream();

      string line;
      string buffer = "";
      while (buffer != "Name")
        {
          getline(f, line);
          istringstream iss(line);
          iss >> buffer; 
        }

      //make shifts list
      string systlabel, dummy;
      float systindex, value, error;
      while (getline(f, line))
        {
          istringstream iss(line);
          iss >> systindex >> systlabel  >> value  >> dummy  >> error; 
          
          shtype sh;
          sh.val = value;
          sh.err = error;
          sh.dataid = itl - opts.labels.begin();
          shlist[systlabel].push_back(sh);
        }
      f.close();
    }

  int s = 0;   //plot number
  //loop on shifts plots
  while (shlist.size() > 0)
    {
      int nshifts = min(opts.spp, (int)shlist.size());
      
      //make canvas
      char numb[15];
      sprintf(numb, "shifts_%d", s);

      //fractional shift heigth
      double shiftheigth = opts.shgth / 1200.;
      //fractional offset heigth
      double offset = 0.4 * shiftheigth;

      int minshifts = (int)((0.085 + 0.025* dirs.size()) / shiftheigth);
  
      //equivalent number of shifts for pad heigth
      //limit the minimum to minshifts, the maximum to the pad heigth divided by shift points
      //      double graphnshifts = max(minshifts, min(nshifts, (int)((1 - tmarg - bmarg - 2 * offset) / shiftheigth)));
      double graphnshifts = max((double)minshifts, min((double)nshifts, (1 - tmarg - bmarg - 2 * offset) / shiftheigth));

      //fractional heigth of pad
      double heigthfact = graphnshifts * shiftheigth + 2 * offset;
      float reduction = max(0., (1 - (tmarg + bmarg + heigthfact)));
      float topreduction = (float)abs(min(0, nshifts - minshifts)) * shiftheigth;
      double poffset = offset * (float)graphnshifts / (float)max(minshifts, nshifts);

      TCanvas * cnv = new TCanvas(numb, "", 0, 0, 2 * opts.resolution, 2 * opts.resolution);
      shiftscnv.push_back(cnv);
      cnv->Divide(2, 1);
      TPad * main = (TPad*)cnv->GetPad(1);
      main->SetLeftMargin(lmarg+0.02);
      main->SetRightMargin(0.01);
      main->SetTopMargin(tmarg);
      main->SetBottomMargin(bmarg + reduction);
      main->cd();

      float x[dirs.size()][nshifts];
      float xerr[dirs.size()][nshifts];
      float y[dirs.size()][nshifts];
      for (int d = 0; d < dirs.size(); d++)
        for (int s = 0; s < nshifts; s++)
          {
            y[d][s] = -100;
            x[d][s] =  0;
            xerr[d][s] = 0;
          }
      //loop on shifts
      map<string, vector<shtype> >::iterator sit = shlist.begin();
      map<string, vector<shtype> >::iterator sitlast = shlist.begin();
      advance(sitlast, nshifts);
      int i = 0;
      for (shlisttype::iterator sit = shlist.begin(); sit != sitlast; sit++, i++)
        {
          //loop on directories
          for (vector<shtype>::iterator dit = (*sit).second.begin(); dit != (*sit).second.end(); dit++)
            {
              x[(*dit).dataid][i] = (*dit).val;
              xerr[(*dit).dataid][i] = (*dit).err;
              y[(*dit).dataid][i] = i + s*opts.spp + 1;
              //displace
              y[(*dit).dataid][i] += -0.5 + ((float)(*dit).dataid + 1) / (float)(dirs.size() + 1);
            }
        }

      //make legend             
      TLegend * leg2 = new TLegend(0.6, 1-tmarg-0.087-dirs.size()*0.025, 1-0.02, 1-tmarg-0.087);
      leg2->SetFillColor(0);
      leg2->SetBorderSize(0);
      leg2->SetTextAlign(12);
      leg2->SetTextFont(62);
      leg2->SetTextSize(txtsize);

      //make TGraphs             
      int nd = 0;
      for (vector<string>::iterator dit = dirs.begin(); dit != dirs.end(); dit++)
        {
          TGraphErrors * gshift = new TGraphErrors(nshifts, x[nd], y[nd], xerr[nd]);
          if (dirs.size() == 1)
            {
              gshift->SetMarkerSize(4.*opts.resolution/1200.);
              gshift->SetMarkerStyle(8);
            }
          else
            {
              gshift->SetMarkerSize(3.*opts.resolution/1200.);
              gshift->SetMarkerStyle(opts.markers[opts.labels[nd]]);
              gshift->SetLineColor(opts.colors[opts.labels[nd]]);
              gshift->SetMarkerColor(opts.colors[opts.labels[nd]]);
            }
          gshift->SetTitle("");
          gshift->GetYaxis()->SetNdivisions(max(minshifts,nshifts)+1);
          gshift->SetMinimum(s*opts.spp+0.1);
          gshift->SetMaximum(s*opts.spp+max(minshifts,nshifts)+0.9);
          gshift->GetXaxis()->SetLabelFont(62);
          gshift->GetXaxis()->SetLabelSize(txtsize);
          gshift->GetYaxis()->SetLabelFont(62);
          gshift->GetYaxis()->SetLabelSize(txtsize);
          gshift->GetXaxis()->SetNdivisions(505);
          gshift->GetXaxis()->Set(100, -3, 7);
          gshift->GetXaxis()->SetTickLength(0.01);
          if (nd == 0)
            gshift->Draw("AP");
          else
            gshift->Draw("P same");
          leg2->AddEntry(gshift, (opts.labels[nd]).c_str(), "p");
          nd++;
        }
 
      TLine *one = new TLine(1, s*opts.spp+0.1, 1, s*opts.spp+nshifts+0.9);
      one->SetLineStyle(2);
      TLine *minusone = new TLine(-1, s*opts.spp+0.1, -1, s*opts.spp+nshifts+0.9);
      minusone->SetLineStyle(2);
      one->Draw();
      minusone->Draw();

           //horizontal grid
          if (dirs.size() > 1)
            for (int h = 0; h <= nshifts; h++)
              {
                TLine *hgrid = new TLine(-3, s*opts.spp+h+0.5, 7, s*opts.spp+h+0.5);
                hgrid->SetLineStyle(3);
                hgrid->SetLineWidth(0.5);
                hgrid->Draw();
            }

      leg2->Draw();

      TPad * legend = (TPad*)cnv->GetPad(2);
      legend->SetLeftMargin(0.01);
      legend->SetRightMargin(rmarg);
      legend->SetTopMargin(tmarg);
      legend->SetBottomMargin(bmarg + reduction);
      legend->cd();
      TPaveText * leg = new TPaveText(0, bmarg+reduction+poffset, 1, 1-tmarg-topreduction-poffset);
      leg->SetFillColor(0);
      leg->SetBorderSize(0);
      leg->SetTextAlign(12);
      leg->SetTextFont(62);
      leg->SetTextSize(txtsize);
      i = nshifts + s*opts.spp;
      for (shlisttype::iterator slab = sitlast; slab != shlist.begin();)
        {
          slab--;
          char num[10];
          sprintf (num, "%d", i);
          string label = (string) num + "  " + slab->first.c_str();
          leg->AddText(label.c_str());
          i--;
        }
      leg->Draw();

      //delete plotted shifts
      shlist.erase(shlist.begin(), sitlast);
      s++;

      cnv->cd();
      if (opts.drawlogo)
	{
	  TPad * xfitterlogo = DrawLogo();
	  float dx = 0.1183 * 1.5;
	  float dy = 0.0744 * 1.5;
	  float xl, yl;
	  xl = 0.5-0.01-0.01;
	  yl = 1-tmarg-0.015;

	  xfitterlogo->SetPad(xl-dx/2, yl-dy/2, xl, yl);
	  xfitterlogo->Draw();
	}
      DrawLabels("ur half");
    }

  return shiftscnv;
}
