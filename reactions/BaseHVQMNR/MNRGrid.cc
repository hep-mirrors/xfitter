#include <MNRGrid.h>
#include "xfitter_cpp.h"
#include <TMath.h>
#include <TSpline.h>

namespace MNR
{
  Grid::Grid()
  {
    fNL   = 0;
    fNY   = 0;
    fNW   = 0;
    fL    = NULL;
    fAs   = NULL;
    fMr   = NULL;
    fY    = NULL;
    fW    = NULL;
    fBW   = NULL;
    fCS   = NULL;
    fNContr = 1;
    fContr = new MNRContribution*[1];
    fContr[0] = new MNRContribution(11111);
    fCS = new double***[fNContr];
    for(int i = 0; i < fNContr; i++) fCS[i] = NULL;
  }

  Grid::Grid(int ncontr, MNRContribution** contr)
  {
    fNL   = 0;
    fNY   = 0;
    fNW   = 0;
    fL    = NULL;
    fAs   = NULL;
    fMr   = NULL;
    fY    = NULL;
    fW    = NULL;
    fBW   = NULL;
    fCS   = NULL;
    fNContr = ncontr;
    fContr = new MNRContribution*[fNContr];
    for(int i = 0; i < fNContr; i++) fContr[i] = new MNRContribution(*contr[i]);
    fCS = new double***[fNContr];
    for(int i = 0; i < fNContr; i++) fCS[i] = NULL;
  }

  Grid::~Grid()
  {
    //printf("OZ Grid::~Grid()\n");
    if(fL) delete fL;
    if(fAs) delete fAs;
    if(fMr) delete fMr;
    if(fY) delete fY;
    if(fW) delete fW;
    if(fBW) delete fBW;
    if(fCS) {
      for(int c = 0; c < fNContr; c++)
      {
        for(int i = 0; i < fNL; i++)
        {
          for(int j = 0; j < fNY; j++)
          {
            delete fCS[c][i][j];
          }
          delete fCS[c][i];
        }
        delete fCS[c];
      }
      delete fCS;
    }
    if(fContr)
    {
      for(int i = 0; i < fNContr; i++)
      {
        delete fContr[i];
      }
      delete fContr;
    }
  }

  void Grid::NonPhys()
  {
    for(int i = 0; i < fNL; i++) this->NonPhys(i);
  }

  void Grid::NonPhys(int bpt)
  {
    for(int i = 0; i < fNY; i++)
      for(int j = 0; j < fNW; j++)
        for(int k = 0; k < fNContr; k++)
          this->CS(k, bpt, i, j) = -1.0; // negative non-physical value
  }

  void Grid::Zero()
  {
    for(int bpt = 0; bpt < fNL; bpt++)
      for(int i = 0; i < fNY; i++)
        for(int j = 0; j < fNW; j++)
          for(int k = 0; k < fNContr; k++)
            this->CS(k, bpt, i, j) = 0;
  }

  void Grid::SetL(int n, double minpt, double maxpt, double xm)
  {
    double power = 0.25;
    fNL = n;
    if(fL) delete fL;
    fL = new double[fNL];
    double xm2 = xm * xm;
    double minpower = TMath::Power(minpt,power);
    double maxpower = TMath::Power(maxpt,power);
    double steppower = (maxpower - minpower) / fNL;
    for(int i = 0; i < fNL; i++)
    {
      double pt = TMath::Power(minpower + i * steppower, 1.0 / power);
      fL[i] = xm2 / (xm2 + pt * pt);
    }
    // array for alpha_s values
    if(fAs) delete fAs;
    fAs = new double[fNL];
    // array for mu_r values
    if(fMr) delete fMr;
    fMr = new double[fNL];
  }

  void Grid::FillPt(double* ptall, double xm)
  {
    double xm2 = xm * xm;
    for(int i = 0; i < fNL; i++) ptall[i] = TMath::Sqrt(xm2 / fL[i] - xm2);
  }

  void Grid::SetY(int n, double min, double max)
  {
    if(!fNL || !fL) {
      std::string str = "F: ERROR in Grid::SetY(): first call Grid::SetL(), then Grid::SetY()\n";
      hf_errlog_(16123010, str.c_str(), str.length());
    }
    if(n < 2 && min != max) {
      char str[256];
      sprintf(str, "F: ERROR in Grid::SetY(): n %d < 2 and min %e != max %e\n", n, min, max);
      hf_errlog_(16123010, str, std::string(str).length());
    }
    fNY = n;
    if(fY) delete fY;
    fY = new double[fNY];
    double step = (max - min) / (n - 1);
    for(int i = 0; i < n; i++) fY[i] = min + step * i;
    for(int c = 0; c < fNContr; c++)
    {
      if(fCS[c]) {
        for(int i = 0; i < fNL; i++) if(fCS[c][i]) delete fCS[c][i];
        delete fCS[c];
      }
      fCS[c] = new double**[fNL];
      for(int i = 0; i < fNL; i++)
      {
        fCS[c][i] = new double*[fNY];
        for(int j = 0; j < fNY; j++) fCS[c][i][j] = NULL;
      }
    }
  }

  void Grid::SetW(int n, double min/* = 0.0*/, double max/* = 500.0*/)
  {
    if(!fNY || !fY)
    {
      std::string str = "F: ERROR in Grid::SetW(): first call Grid::SetY(), then Grid::SetW()\n";
      hf_errlog_(16123010, str.c_str(), str.length());
    }
    fNW = n;
    if(fW) delete fW;
    fW = new double[fNW];
    if(fBW) delete fBW;
    fBW = new double[fNW+1];
    double step = (max - min) / fNW;
    for(int i = 0; i <= fNW; i++)
    {
      fBW[i] = min + step * i * i / fNW;
      if(i == 0) continue;
      fW[i-1] = (fBW[i] + fBW[i-1]) / 2;
    }
    for(int c = 0; c < fNContr; c++)
      for(int i = 0; i < fNL; i++)
        for(int j = 0; j < fNY; j++)
        {
          if(fCS[c][i][j]) delete fCS[c][i][j];
          fCS[c][i][j] = new double[fNW];
        }
  }

  void Grid::SetW(double b1, double b2)
  {
    if(!fNY || !fY)
    {
      std::string str = "F: ERROR in Grid::SetW(): first call Grid::SetY(), then Grid::SetW()\n";
      hf_errlog_(16123010, str.c_str(), str.length());
    }
    fNW = 3;
    if(fW) delete fW;
    fW = new double[fNW];
    fW[0] = b1 / 2.;
    fW[1] = (b2 - b1) / 2.;
    fW[2] = b2 * 2.;
    if(fBW) delete fBW;
    fBW = new double[fNW+1];
    fBW[0] = 0.0;
    fBW[1] = b1;
    fBW[2] = b2;
    fBW[3] = b2 * 3.;
    for(int c = 0; c < fNContr; c++)
      for(int i = 0; i < fNL; i++)
        for(int j = 0; j < fNY; j++)
        {
          if(fCS[c][i][j]) delete fCS[c][i][j];
          fCS[c][i][j] = new double[fNW];
        }
  }

  int Grid::FindWBin(double w)
  {
    for(int i = 0; i < fNW; i++)
        if(w < fBW[i+1] && w > fBW[i]) return i;
    return fNW - 1;
  }

  void Grid::Print(double xm)
  {
    double xm2 = xm * xm;
    for(int c = 0; c < fNContr; c++)
    {
      for(int bpt = 0; bpt < fNL; bpt++)
      {
        double mt2 = xm2 / fL[bpt];
        double pt = TMath::Sqrt(mt2 - xm2);
        for(int by = 0; by < fNY; by++)
          for(int bw = 0; bw < fNW; bw++)
            printf("pt: %e  y: %e  w: %e  CS: %e\n", pt, fY[by], fW[bw], fCS[c][bpt][by][bw]);
      }
    }
  }

  void Grid::InterpolateGrid(Grid* gridorig, Grid* gridtrg, double mq)
  {
    // Get pT array of target grid
    double pt[gridtrg->fNL];
    gridtrg->FillPt(pt, mq);
    // Transform this pT array into array of L = m^2 / (m^2 + pT^2)
    double mq2 = mq * mq;
    int nltrg = gridtrg->NL();
    double* ltrg = gridtrg->LPtr();
    for(int i = 0; i < nltrg; i++) ltrg[i] = mq2 / (mq2 + pt[i] * pt[i]);
    // For spline: prepare L array of original grid in reversed order
    int nlorig = gridorig->NL();
    double* lorig = gridorig->LPtr();
    double spline_x[nlorig], spline_y[nlorig];
    for(int i = 0; i < nlorig; i++)
      spline_x[nlorig-1-i] = lorig[i];
    // Loop over contributions
    for(int c = 0; c < gridorig->GetNContr(); c++)
      // Loop over y bins
      for(int y = 0; y < gridorig->NY(); y++)
        // Loop over W bins
        for(int w = 0; w < gridorig->NW(); w++)
        {
          // For spline: prepare X-section array of original grid in reversed order
          for(int l = 0; l < nlorig; l++) spline_y[nlorig-1-l] = gridorig->CS(c,l,y,w);
          // Spline interpolation
          TSpline3 spline("", spline_x, spline_y, nlorig);
          for(int l = 0; l < nltrg; l++) gridtrg->CS(c,l,y,w) = spline.Eval(ltrg[l]);
        }
    // interpolate as and mu_r
    double spline_y_as[nlorig], spline_y_mr[nlorig];
    for(int l = 0; l < nlorig; l++)
    {
      spline_y_as[nlorig-1-l] = gridorig->AlphaS(l);
      spline_y_mr[nlorig-1-l] = gridorig->MuR(l);
    }
    TSpline3 spline_as("", spline_x, spline_y_as, nlorig);
    TSpline3 spline_mr("", spline_x, spline_y_mr, nlorig);
    for(int l = 0; l < nltrg; l++)
    {
      gridtrg->AlphaS(l) = spline_as.Eval(ltrg[l]);
      gridtrg->MuR(l) = spline_mr.Eval(ltrg[l]);
    }
  }

  // Different mass schemes:
  // flag = 0: MSbar mass scheme with mu = mu_R
  // flag = 1: MSbar mass scheme with mu = mu_R and l != 0
  // flag = 2: MSR mass scheme with R provided (nf should be provided - number of flavours)
  void Grid::InterpolateGrid(Grid *gridorig, Grid *gridtrg, double mq, Grid *gridorig_LO_massUp, double mq_masUp, Grid *gridorig_LO_massDown, double mq_masDown, int flag, double* R, int* nf, double mq_mu0)
  {
    double mqDiff = mq_masUp - mq_masDown;
    // Get pT array of target grid
    double pt[gridtrg->fNL];
    gridtrg->FillPt(pt, mq);
    // Transform this pT array into array of L = m^2 / (m^2 + pT^2)
    double mq2 = mq * mq;
    double mq_masUp2 = mq_masUp * mq_masUp;
    double mq_masDown2 = mq_masDown * mq_masDown;
    int nltrg = gridtrg->NL();
    double* ltrg = gridtrg->LPtr();
    for(int i = 0; i < nltrg; i++) ltrg[i] = mq2 / (mq2 + pt[i] * pt[i]);
    // For spline: prepare L array of original grid in reversed order
    int nlorig = gridorig->NL();
    double* lorig = gridorig->LPtr();
    double* lorig_LO_massUp = gridorig_LO_massUp->LPtr();
    double* lorig_LO_massDown = gridorig_LO_massDown->LPtr();
    double spline_x[3][nlorig], spline_y[5][nlorig];
    for(int i = 0; i < nlorig; i++)
    {
      spline_x[0][i] = mq2 / lorig[i] - mq2;
      spline_x[1][i] = mq_masUp2 / lorig_LO_massUp[i] - mq_masUp2;
      spline_x[2][i] = mq_masDown2 / lorig_LO_massDown[i] - mq_masDown2;
    }
    // Loop over contributions
    for(int c = 0; c < gridorig->GetNContr(); c++)
      // Loop over y bins
      for(int y = 0; y < gridorig->NY(); y++)
        // Loop over W bins
        for(int w = 0; w < gridorig->NW(); w++)
        {
          // For spline: prepare X-section array of original grid in reversed order
          for(int l = 0; l < nlorig; l++)
          {
            spline_y[0][l] = gridorig->CS(c,l,y,w);
            spline_y[1][l] = gridorig_LO_massUp->CS(c,l,y,w);
            spline_y[2][l] = gridorig_LO_massDown->CS(c,l,y,w);
            spline_y[3][l] = gridorig->AlphaS(l);
            spline_y[4][l] = gridorig->MuR(l);
          }
          // Spline interpolation
          TSpline3 spline("", spline_x[0], spline_y[0], nlorig);
          TSpline3 spline_LO_massUp("", spline_x[1], spline_y[1], nlorig);
          TSpline3 spline_LO_massDown("", spline_x[2], spline_y[2], nlorig);
          TSpline3 spline_as("", spline_x[0], spline_y[3], nlorig);
          TSpline3 spline_mr("", spline_x[0], spline_y[4], nlorig);
          for(int l = 0; l < nltrg; l++)
          {
            double pt2 = pt[l] * pt[l];
            double xsecOld = spline.Eval(pt2);
            double xsecOld_LO_massUp = spline_LO_massUp.Eval(pt2);
            double xsecOld_LO_massDown = spline_LO_massDown.Eval(pt2);
            double delta = 0.0;
            if(flag == 0 || flag == 1 || flag == 11)
            {
              double as = spline_as.Eval(pt2);
              double mr = spline_mr.Eval(pt2);
              double d1 = 4.0 / 3.0;
              if(flag == 1)
                d1 += 2.0 * TMath::Log(mr / mq);
              else if(flag == 11)
                //d1 += 2.0 * TMath::Log(2.*mr / mq);
                d1 += 2.0 * TMath::Log(mq_mu0 / mq);
              delta = as / TMath::Pi() * d1 * mq * (xsecOld_LO_massUp - xsecOld_LO_massDown) / mqDiff;
            }
            else if(flag == 2)
            {
              double as = spline_as.Eval(*R);
              double b0 = 11.0 - 2.0 / 3.0 * (*nf);
              double a1 = 2 * b0 * 0.348;
              delta = (*R) * a1 * as / (4.0 * TMath::Pi()) * (xsecOld_LO_massUp - xsecOld_LO_massDown) / mqDiff;
            }
            double xsecNew = xsecOld + delta;
            gridtrg->CS(c,l,y,w) = xsecNew;
          }
        }
  }

  void Grid::TransformGridToMSbarMassScheme(Grid *grid, Grid *gridLOMassUp, Grid *gridLOMassDown, double mq, double mqDiff)
  {
    // Loop over contributions
    for(int c = 0; c < grid->GetNContr(); c++)
    {
      // Loop over y bins
      for(int y = 0; y < grid->NY(); y++)
      {
        // Loop over W bins
        for(int w = 0; w < grid->NW(); w++)
        {
          // Loop over pT (L) bins
          for(int l = 0; l < grid->NL(); l++)
          {
            double as = grid->AlphaS(l);
            //printf("as = %f\n", as);
            double d1 = 4.0 / 3.0;
            double xsecNew = grid->CS(c,l,y,w) + as / TMath::Pi() * d1 * mq * (gridLOMassUp->CS(c,l,y,w) - gridLOMassDown->CS(c,l,y,w)) / (2 * mqDiff);
            if(y == 20)
              printf("c,y,w,l,g,u,d: %d %d %d %d %f[%f] %f[%f] %f[%f]  -->  %f [%.2f]\n", c, y, w, l, grid->CS(c,l,y,w), grid->LPtr()[l], gridLOMassUp->CS(c,l,y,w), gridLOMassUp->LPtr()[l], gridLOMassDown->CS(c,l,y,w), gridLOMassDown->LPtr()[l], xsecNew, xsecNew / grid->CS(c,l,y,w) * 100);
            grid->CS(c,l,y,w) = xsecNew;
          }
        }
      }
    }
  }
}
