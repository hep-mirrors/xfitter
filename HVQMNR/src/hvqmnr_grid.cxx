#include <fortran_wrapper.h>
#include <hvqmnr_grid.h>
#include <TMath.h>
#include <TSpline.h>


namespace HVQMNR
{
  Grid::Grid() 
  {
    fNL   = 0;
    fNY   = 0;
    fNW   = 0;
    fL    = NULL;
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
    if(fL) delete fL;
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
  }

  void Grid::FillPt(double* ptall, double xm) 
  {
    double xm2 = xm * xm;
    for(int i = 0; i < fNL; i++) ptall[i] = TMath::Sqrt(xm2 / fL[i] - xm2);
  }

  void Grid::SetY(int n, double min, double max)
  {
    if(!fNL || !fL) {
      printf("ERROR in Grid::SetY(): first call Grid::SetL(), then Grid::SetY()\n");
      hf_stop_();
    }
    if(n < 2 && min != max) {
      printf("ERROR in Grid::SetY(): n %d < 2 and min %e != max %e\n", n, min, max);
      hf_stop_();
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
      printf("ERROR in Grid::SetW(): first call Grid::SetY(), then Grid::SetW()\n");
      hf_stop_();
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
      printf("ERROR in Grid::SetW(): first call Grid::SetY(), then Grid::SetW()\n");
      hf_stop_();
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
  }
}
