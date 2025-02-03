#include <MNRGrid.h>
#include <MNRFrag.h>
#include "xfitter_cpp.h"
#include <TMath.h>
#include <TF2.h>
#include <TH2D.h>
#include "Math/IntegratorOptions.h"

// Interface to xFitter FORTRAN routines
extern "C"
{
  void hf_stop_();
}

namespace MNR
{
  Frag::Frag()
  {
    //printf("OZ Frag::Frag()\n");
    fNout  = 0;
    fBnz   = 0;
    fCGnpt = 0;
    fCGny  = 0;
    fCGpt  = NULL;
    fCGy   = NULL;
    fCGptf = NULL;
    fCGyf  = NULL;
    fZc    = NULL;
    fZw    = NULL;
    bDebug = false;
    bFirst = true;
    fLastxm = -1.0;
    fNRecalc = 0;
  }

  Frag::~Frag()
  {
    //printf("OZ Frag::~Frag()\n");
    this->ClearZ();
    this->ClearPrecalc();
    for(unsigned int f = 0; f < fFFrag.size(); f++) if(fFFrag[f]) delete fFFrag[f];
  }

  void Frag::ClearZ()
  {
    if(fBnz)
    {
      delete fZc;
      fZc = NULL;
      delete fZw;
      fZw = NULL;
      for(int f = 0; f < fNout; f++)
      {
        delete fWz[f];
        fWz[f] = NULL;
      }
    }
    fBnz = 0;
  }

  void Frag::ClearPrecalc()
  {
    if(fCGnpt && fCGny)
    {
      delete fCGpt;
      fCGpt = NULL;
      delete fCGy;
      fCGy = NULL;
      for(int bpt = 0; bpt < fCGnpt; bpt++)
      {
        delete fCGptf[bpt];
        fCGptf[bpt] = NULL;
        for(int by = 0; by < fCGny; by++)
        {
          for(int bz = 0; bz < fBnz; bz++)
          {
            delete fCGyf[bpt][by][bz];
            fCGyf[bpt][by][bz] = NULL;
          }
          delete fCGyf[bpt][by];
          fCGyf[bpt][by] = NULL;
        }
        delete fCGyf[bpt];
        fCGyf[bpt] = NULL;
      }
      delete fCGptf;
      fCGptf = NULL;
      delete fCGyf;
      fCGyf = NULL;
    }
    fCGnpt = fCGny = 0;
  }

  void Frag::SetNz(int nz)
  {
    if(nz < 1)
    {
      char str[256];
      sprintf(str, "F: ERROR in Frag::SetNz(): nz %d < 1\n", nz);
      hf_errlog_(16123010, str, std::string(str).length());
    }
    this->ClearZ();
    fBnz = nz;
    fZc = new double[fBnz];
    fZw = new double[fBnz];
    double zprev, znext;
    zprev = 0.0;
    for(int bz = 0; bz < fBnz+1; bz++)
    {
      znext = 1. * bz / fBnz;
      znext = TMath::Power(znext, 0.75);
      if(bz == 0) continue;
      fZc[bz-1] = (znext + zprev) / 2;
      fZw[bz-1] = znext - zprev;
      zprev = znext;
    }
  }

  void Frag::AddOut(TF1* ffrag, double M)
  {
    if(!fBnz)
    {
      std::string str = "F: ERROR in Frag::AddOut(): first call Frag::SetNz()\n";
      hf_errlog_(16123010, str.c_str(), str.length());
    }
    fFFrag.push_back(ffrag);
    // If negative mass is provided, use heavy-quark mass instead
    if(M < 0) fMh2.push_back(-1.0);
    else fMh2.push_back(M * M);
    fWz.push_back(new double[fBnz]);
    fNout++;
  }

  void Frag::Precalc(Grid* grid, double xm)
  {
    fNRecalc++;
    // Fast variables
    this->ClearPrecalc();
    fCGnpt = grid->NL();
    double* p_l = grid->LPtr();
    fCGny = grid->NY();
    double* p_y = grid->YPtr();
    fCGpt = new double[fCGnpt];
    fCGy = new double[fCGny];
    fCGptf = new double*[fCGnpt];
    fCGyf = new double***[fCGnpt];
    double xm2 = xm * xm;
    double M2[fNout];
    for(int f = 0; f < fNout; f++)
    {
      M2[f] = fMh2[f];
      if(fMh2[f] < 0.0) M2[f]=xm2;
    }
    double shy[fCGny];
    for(int by = 0; by < fCGny; by++)
    {
      fCGy[by] = p_y[by];
      shy[by] = TMath::SinH(fCGy[by]);
    }
    // Loop over pT (internally L)
    for(int bpt = 0; bpt < fCGnpt; bpt++)
    {
      double l2 = p_l[bpt];
      double mt2 = xm2 / l2;
      double mt = TMath::Sqrt(mt2);
      double pt2 = mt2 - xm2;
      fCGpt[bpt] = TMath::Sqrt(pt2);
      if(mt2 < xm2) continue;
      fCGptf[bpt] = new double[fBnz];
      fCGyf[bpt] = new double**[fCGny];
      // Loop over y
      for(int by = 0; by < fCGny; by++)
      {
        double pl = mt * shy[by];
        fCGyf[bpt][by] = new double*[fBnz];
        // Loop over z
        for(int bz = 0; bz < fBnz; bz++)
        {
          if(by == 0) fCGptf[bpt][bz] = fCGpt[bpt] * fZc[bz];
          fCGyf[bpt][by][bz] = new double[fNout];
          double plf = pl * fZc[bz];
          // Loop over final states
          for(int f = 0; f < fNout; f++)
          {
            double Mt = TMath::Sqrt(M2[f] + fCGptf[bpt][bz] * fCGptf[bpt][bz]);
            fCGyf[bpt][by][bz][f] = TMath::ASinH(plf / Mt);
          }
        }
      }
    }
  }

  void Frag::CalcCS(Grid* grid, double xm, std::vector<TH2D*> hcs)
  {
    // If it is first run or heavy-quark mass changed, recalculation is needed
    if(bFirst || fLastxm != xm)
    {
      this->Precalc(grid, xm);
      fLastxm = xm;
    }
    if(bFirst) bFirst=0;

    // Prepare output histograms and fast variables
    int ncontr = grid->GetNContr();
    MNRContribution* contr[ncontr];
    int nw = grid->NW();
    double* p_w = grid->WPtr();
    double wgrid[ncontr][fCGnpt][fCGny][nw];
    double* harray[ncontr][fNout];
    for(int c = 0; c < ncontr; c++)
    {
      contr[c] = grid->GetContr(c);
      for(int f = 0; f < fNout; f++)
      {
        TH2D* h = hcs[f * ncontr + c];
        h->Reset();
        harray[c][f] = h->GetArray();
      }
      for(int bpt = 0; bpt < fCGnpt; bpt++)
        for(int by = 0; by < fCGny; by++)
          for(int bw = 0; bw < nw; bw++)
            wgrid[c][bpt][by][bw] = grid->CS(c, bpt, by, bw);
    }

    // Check heavy-quark mass for nan
    if(xm != xm) return;
    double wfz[fNout][fBnz][nw];
    for(int f = 0; f < fNout; f++)
    {
      // Parton level
      if(fFFrag[f] == 0) continue;
      // 1D fragmentation function
      if(fFFrag[f]->ClassName() == TString("TF1"))
      {
        double norm = 0.0;
        for(int bz = 0; bz < fBnz; bz++)
        {
          double z = fZc[bz];
          double wz = fZw[bz] * fFFrag[f]->Eval(z);
          norm += wz;
          for(int bw = 0; bw < nw; bw++) wfz[f][bz][bw]=wz;
        }
        for(int bz = 0; bz < fBnz; bz++)
          for(int bw = 0; bw < nw; bw++)
            wfz[f][bz][bw] /= norm;
      }
      // 2D fragmentation function
      else if(fFFrag[f]->ClassName() == TString("TF2"))
      {
        for(int bw = 0; bw < nw; bw++)
        {
          double w = p_w[bw];
          double norm = 0.0;
          for(int bz = 0; bz < fBnz; bz++)
          {
            double z = fZc[bz];
            double zw = fZw[bz];
            wfz[f][bz][bw] = zw * fFFrag[f]->Eval(z,w);
            norm += wfz[f][bz][bw];
          }
          for(int bz = 0; bz < fBnz; bz++) wfz[f][bz][bw] /= norm;
        }
      }
      else
      {
        char str[256];
        sprintf(str, "F: ERROR in Frag::CalcCS(): ff[%d] does not belong to TF1 or TF2\n", f);
        hf_errlog_(16123010, str, std::string(str).length());
      }
    }

    // Loop over pT
    for(int bpt = 0; bpt < fCGnpt - 1; bpt++)
    {
      int bptn = bpt + 1;
      if(bDebug) if(bpt % 10 == 0) printf("Frag::CalcCS() 1st dim: %3d from %3d\n", bpt, fCGnpt - 1);
      double pt1 = fCGpt[bpt];
      if(pt1 != pt1) continue;
      double pt2 = fCGpt[bptn];
      double pt_w = pt2 - pt1;
      // Loop over y
      for(int by = 0; by < fCGny - 1; by++)
      {
        int byn = by + 1;
        double y1 = fCGy[by];
        double y2 = fCGy[byn];
        double y_w = y2 - y1;
        double pty_w = pt_w * y_w;
        // Loop over z
        for(int bz = 0; bz < fBnz; bz++)
        {
          // Loop over final states
          for(int f = 0; f < fNout; f++)
          {
            double wf, pth1, pth2, yh11, yh12, yh21, yh22;
            if(!fFFrag[f])
            { // parton level: no rescaling
              if(bz != 0) continue;
              wf = 1.0;
              pth1 = pt1;
              pth2 = pt2;
              yh11 = yh21 = y1;
              yh12 = yh22 = y2;
            }
            else
            {
              wf = 0.0;
              pth1 = fCGptf[bpt][bz];
              pth2 = fCGptf[bptn][bz];
              yh11 = fCGyf[bpt][by][bz][f];
              yh12 = fCGyf[bpt][byn][bz][f];
              yh21 = fCGyf[bptn][by][bz][f];
              yh22 = fCGyf[bptn][byn][bz][f];
            }
            // Loop over contributions
            for(int c = 0; c < ncontr; c++)
            {
              if(contr[c]->fActive == 0) continue;
              int offset = f * ncontr + c;
              TH2D* h = hcs[offset];
              TAxis* ax = h->GetXaxis();
              int baxnb = ax->GetNbins();
              int bax1 = ax->FindBin(pth1);
              int bax2 = ax->FindBin(pth2);
              if(bax2 == 0 || bax1 == baxnb + 1) continue;
              TAxis* ay = h->GetYaxis();
              int baynb = ay->GetNbins();
              double xw = pth2 - pth1;
              // Loop over hadron pT bins
              for(int bax = bax1; bax <= bax2; bax++)
              {
                double xl = pth1;
                double xh = pth2;
                if(bax > bax1) xl = ax->GetBinLowEdge(bax);
                if(bax < bax2) xh = ax->GetBinUpEdge(bax);
                double xc = (xl + xh) / 2;
                double xh_m_xl_ov_xw = (xh - xl) / xw;
                double xwcoef = xh_m_xl_ov_xw / xw;
                double dx = (xc - pth1) * xwcoef;
                double dxr = (pth2 - xc) * xwcoef;
                double yxwcoef = xwcoef * xw;
                double yh1 = (yh11 * dxr + yh21 * dx) / yxwcoef;
                double yh2 = (yh12 * dxr + yh22 * dx) / yxwcoef;
                double yw = yh2 - yh1;
                int bay1 = ay->FindBin(yh1);
                int bay2 = ay->FindBin(yh2);
                if(bay2 == 0 || bay1 == baynb + 1) continue;
                // Loop over hadron y bins
                for(int bay = bay1; bay <= bay2; bay++)
                {
                  double yl = yh1;
                  double yh = yh2;
                  if(bay>bay1) yl = ay->GetBinLowEdge(bay);
                  if(bay<bay2) yh = ay->GetBinUpEdge(bay);
                  double yc = (yl + yh) / 2;
                  double ywcoef = (yh - yl) / (yw * yw);
                  double dy = (yc - yh1) * ywcoef;
                  double dyr = (yh2 - yc) * ywcoef;
                  int binxy = h->GetBin(bax, bay);
                  // Loop over W
                  for(int bw = 0; bw < nw; bw++)
                  {
                    if(fFFrag[f]) wf = wfz[f][bz][bw]; // not parton level
                    // Smooth final contrbution as linear 2D function
                    double w_parton = pty_w * (wgrid[c][bpt][by][bw] * dxr * dyr +
                                               wgrid[c][bpt][byn][bw] * dxr * dy +
                                               wgrid[c][bptn][by][bw] * dx * dyr +
                                               wgrid[c][bptn][byn][bw] * dx * dy);
                    harray[c][f][binxy] += w_parton * wf;
                  } // W
                } // hadron y
              } // hadron pT
            } // contributions
          } // final states
        } // z
      } // y
    } // pT
  }

  TF1* Frag::GetFragFunction(int f, const char* meson, double par, double* mean/* = 0*/)
  {
    TF1* f_meson = NULL;
    int f1 = f / 10;
    int f2 = f % 10;
    if(f1 == 0)
    {
      if(f2 == 0)
      { // Kartvelishvili
        if(TString(meson) == TString("dzero")) f_meson = new TF1("f_kar_dzero", Frag::kar_dzero, 0., 1., 2);
        else if(TString(meson) == TString("dch")) f_meson = new TF1("f_kar_dch", Frag::kar_dch, 0., 1., 2);
        else f_meson = new TF1("f_kar", Frag::kar, 0., 1., 2);
      }
      else if(f2==1)
      { // BCFY
        if(TString(meson) == TString("dzero")) f_meson = new TF1("f_bcfy_dzero", Frag::bcfy_dzero, 0., 1., 2);
        else if(TString(meson)==TString("dch")) f_meson = new TF1("f_bcfy_dch", Frag::bcfy_dch, 0., 1., 2);
        else f_meson = new TF1("f_bcfy", Frag::bcfy_v, 0., 1., 2);
      }
      else if(f2==2)
      { // Peterson
        f_meson = new TF1("f_pet", Frag::pet, 0., 1., 2);
      }
    }
    else if(f1 == 1)
    {
      if(f2 == 0)
      { // Kartvelishvili "Misha-style"
        if(TString(meson) == TString("dzero")) f_meson = new TF2("f_karw_dzero", Frag::karw_dzero, 0., 1., 0., 10000., 3);
        else if(TString(meson) == TString("dch")) f_meson = new TF2("f_karw_dch", Frag::karw_dch, 0., 1., 0., 10000., 3);
        else f_meson = new TF2("f_karw", Frag::karw, 0., 1., 0., 10000., 3);
      }
    }
    else if(f1 == 2)
    {
      if(f2 == 0)
      { // Kartvelishvili step
        if(TString(meson) == TString("dzero")) f_meson = new TF2("f_karstep_dzero", Frag::karstep_dzero, 0., 1., 0., 10000., 3);
        else if(TString(meson) == TString("dch")) f_meson = new TF2("f_karstep_dch", Frag::karstep_dch, 0., 1., 0., 10000., 3);
        else f_meson = new TF2("f_karstep", Frag::karstep, 0., 1., 0., 10000., 3);
      }
    }
    if(!f_meson)
    {
      char str[256];
      sprintf(str, "F: ERROR in Frag::GetFragFunction(): unknown f %d\n", f);
      hf_errlog_(16123010, str, std::string(str).length());
    }
    f_meson->SetNpx(10000);
    for(int p = 0; p < f_meson->GetNpar(); p++) f_meson->SetParameter(p, 1.0);
    f_meson->SetParameter(1, par);
    // For 1D fragmentation function determine mean, if needed
    if(f_meson->ClassName() == TString("TF1"))
    {
      double integral = f_meson->Integral(0.,1.);
      if (integral == 0.) // SZ 2024.07.01 this happens with ROOT Version: 6.32.02 on naf el9: need to change integration method
      {
        std::string default_mathod = ROOT::Math::IntegratorOneDimOptions::DefaultIntegrator();
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("GAUSS");
        integral = f_meson->Integral(0.,1.);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator(default_mathod.c_str());
      }
      f_meson->SetParameter(0, 1. / integral);
      if(mean) *mean = f_meson->Mean(0.0, 1.0);
    }
    return f_meson;
  }

  double Frag::bcfy_v(double* x, double* p)
  {
    return 3.0 * p[0] * (p[1] * x[0] * TMath::Power((1. - x[0]), 2.) * TMath::Power((1. - (1. - p[1]) * x[0]), -6.) *
           (2. - 2. * (3. - 2 * p[1]) * x[0] + 3. * (3. - 2. * p[1] + 4. * p[1] * p[1]) * x[0] * x[0] - 2. * (1. - p[1]) *
           (4. - p[1] + 2 * p[1] * p[1]) * x[0] * x[0] * x[0] + (1. - p[1]) * (1. - p[1]) * (3. - 2. * p[1] +
           2. * p[1] * p[1]) * x[0] * x[0] * x[0] * x[0]));
  }

  double Frag::bcfy_v_prim(double* x, double* p)
  {
    if(x[0] > (fM_dzero / fM_dstar)) return 0.;
    double newx = fM_dstar /fM_dzero * x[0];
    return (fM_dstar / fM_dzero) * bcfy_v(&newx, p);
  }

  double Frag::bcfy_p(double* x, double* p)
  {
    return p[0] * (p[1] *x[0] * TMath::Power((1. - x[0]), 2.) * TMath::Power((1. - (1. - p[1]) * x[0]), -6.) *
           (6. - 18. * (1. - 2 * p[1]) * x[0] + (21. - 74. * p[1] + 68. * p[1] * p[1]) * x[0] * x[0] - 2. *
           (1. - p[1]) * (6. - 19. * p[1] + 18. * p[1] *p[1]) * x[0] * x[0] * x[0] + 3. * (1. - p[1]) * (1. - p[1]) *
           (1. - 2. * p[1] + 2. * p[1] * p[1]) * x[0] * x[0] * x[0] * x[0]));
  }

  double Frag::bcfy_dzero(double* x, double* p)
  {
    return p[0] * (0.168 * bcfy_p(x, p) + 0.390 * bcfy_v_prim(x, p));
  }

  double Frag::kar_dzero(double* x, double* p)
  {
    return p[0] * (0.168 * kar(x, p) + 0.390 * kar_prim(x, p));
  }

  double Frag::karw_dzero(double* x, double* p)
  {
    return p[0] * (0.168 * karw(x, p) + 0.390 * karw_prim(x, p));
  }

  double Frag::karstep_dzero(double* x, double* p)
  {
    return p[0] * (0.168 * karstep(x, p) + 0.390 * karstep_prim(x, p));
  }

  double Frag::bcfy_dch(double* x, double* p)
  {
    return p[0] * (0.162 * bcfy_p(x, p) + 0.07153 * bcfy_v_prim(x, p));
  }

  double Frag::kar_dch(double* x, double* p)
  {
    return p[0] * (0.162 * kar(x, p) + 0.07153 * kar_prim(x, p));
  }

  double Frag::karw_dch(double* x, double* p)
  {
    return p[0] * (0.162 * karw(x, p) + 0.07153 * karw_prim(x,p));
  }

  double Frag::karstep_dch(double* x, double* p)
  {
    return p[0] * (0.162 * karstep(x, p) + 0.07153 * karstep_prim(x, p));
  }

  double Frag::kar(double* x, double* p)
  {
    return p[0] * TMath::Power(x[0], p[1]) * (1 - x[0]);
  }

  double Frag::kar_prim(double* x, double* p)
  {
    if(x[0] > (fM_dzero / fM_dstar)) return 0.;
    double newx = fM_dstar / fM_dzero * x[0];
    return (fM_dstar / fM_dzero) * kar(&newx, p);
  }

  double Frag::karw(double* x, double* p)
  {
    double alpha = p[1] + p[2] / x[1];
    // prevent very hard form (may lead to numerical problems)
    if(alpha > 100.0) alpha = 100.0;
    // prevent very soft form (may lead to numerical problems)
    if(alpha < 0.1) alpha = 0.1;
    double newp[2] = { p[0], alpha };
    return kar(x, newp);
  }

  double Frag::karw_prim(double* x, double* p)
  {
    if(x[0] > (fM_dzero / fM_dstar)) return 0.;
    double newx[2] = { fM_dstar / fM_dzero * x[0], x[1] };
    return (fM_dstar / fM_dzero) * karw(newx, p);
  }

  double Frag::karstep(double* x, double* p)
  {
    // Values from arXiv:1211.1182
    double alphastep[3][3] = { { 2.36, 2.9, 5.2 }, { 2.67, 3.3, 6.1}, {2.98, 3.7, 7.0} };
    double bin1[3] = { 21.0, 61.0, 101.0 };
    double bin2 = 315.0;
    int var = TMath::Floor(p[1]);
    int bin = TMath::Floor(p[2]);
    double alpha = -1.0;
    if(x[1] < bin1[bin]) alpha = alphastep[var][2]; // bin1
    else if(x[1] < bin2) alpha = alphastep[var][1]; // bin2
    else                 alpha = alphastep[var][0]; // bin3
    double newp[2] = { x[0], alpha };
    return kar(x, newp);
  }

  double Frag::karstep_prim(double* x, double* p)
  {
    if(x[0] > (fM_dzero / fM_dstar)) return 0.;
    double newx = fM_dstar / fM_dzero * x[0];
    return (fM_dstar / fM_dzero) * karstep(&newx, p);
  }

  double Frag::pet(double* x, double* p)
  {
    return p[0] * TMath::Power(x[0], -1.) * TMath::Power(1. - 1./x[0] - p[1]/(1.-x[0]), -2.);
  }

  double Frag::GetHadronMass(const char* meson)
  {
    if(std::string(meson) == "dzero")
      return fM_dzero;
    else if(std::string(meson) == "dch")
      return fM_dch;
    else if(std::string(meson) == "dstar")
      return fM_dstar;
    else if(std::string(meson) == "ds")
      return fM_ds;
    else if(std::string(meson) == "lambdac")
      return fM_lambdac;
    else if(std::string(meson) == "bzero")
      return fM_bzero;
    else if(std::string(meson) == "bch")
      return fM_bch;
    else if(std::string(meson) == "bs")
      return fM_bs;
    else
      return -1.0;
  }


  // Values from PDG
  const double Frag::fM_dzero   = 1.865;
  const double Frag::fM_dch     = 1.867;
  const double Frag::fM_dstar   = 2.010;
  const double Frag::fM_ds      = 1.969;
  const double Frag::fM_lambdac = 2.286;
  const double Frag::fM_bzero   = 5.279;
  const double Frag::fM_bch     = 5.280;
  const double Frag::fM_bs      = 5.367;
}
