#ifndef HVQMNR_GRID_H
#define HVQMNR_GRID_H


namespace HVQMNR
{
  // Structure to store needed contribution for MNR calculation (LO, NLO and gg, qq, qg)
  struct MNRContribution {
    // Constructor
    MNRContribution(int contr)
    {
      fActive = true;
      fLO  = (contr / 10000) % 10;
      fNLO = (contr / 1000)  % 10;
      fgg  = (contr / 100)   % 10;
      fqq  = (contr / 10)    % 10;
      fqg  = (contr / 1)     % 10;
    }

    // LO contribution
    bool fLO;
    // NLO contribution
    bool fNLO;
    // Contribution from gg process
    bool fgg;
    // Contribution from qqbar process
    bool fqq;
    // Contribution from gq process
    bool fqg;
    // Flag to enable or disable all contributions
    // (if fActive = false, all calculations are just omitted)
    bool fActive;
  };


  //
  // Author: Oleksandr Zenaiev (oleksandr.zenaiev@desy.de)
  //
  // Class for grid with binning in pT, y and W (squared energy in parton-parton rest frame) used by MNR and Frag classes
  // (internally binning in L = xm^2 / (xm^2 + pT^2) is used instead of pT, where xm is the heavy-quark mass).
  class Grid {
  // Public methods
  public:
    // Default constructor (one full contrbution: LO+NLO gg+q+qg)
    Grid();
    
    // Construct with specified array of contributions
    Grid(int ncontr, MNRContribution** contr);
    
    // Destructor
    ~Grid();
    
    // Set pT binning (n bins from minpt to maxpt, xm is heavy-quark mass)
    // (internally binning in L = xm^2 / (xm^2 + pT^2))
    void SetL(int n, double minpt, double maxpt, double xm);
    
    // Set y binning (n bins from min to max)
    void SetY(int n, double min, double max);
    
    // Set W (squared energy in parton-parton rest frame), default binning
    // WARNING: it is known feature that the fit might not converge if n != 1 (i.e. if fragmentation depends 
    //          on energy in parton-parton rest frame), especially if the heavy-quark mass is released, 
    //          use this with great caution!
    void SetW(int n = 1, double min = 0.0, double max = 500.0);
    // Set W (squared energy in parton-parton rest frame) three bins, separated by b1 and b2 values
    // (corresponds to the fragmentation set-up used for charm in arXiv:1211.1182)
    // WARNING: it is known feature that the fit might not converge if this option is used, 
    //          especially if the heavy-quark mass is released, use this with great caution!
    void SetW(double b1, double b2);
    
    // Get cross section (by reference) in specified bin
    inline double& CS(int contr, int bl, int by, int bw=0) { return fCS[contr][bl][by][bw]; };

    // Get number of pT (L) bins
    inline int NL() { return fNL; };

    // Get number of y bins
    inline int NY() { return fNY; };

    // Get number of W bins
    inline int NW() { return fNW; };
        
    // Get L binning
    inline double* LPtr() { return fL; };

    // Get y binning
    inline double* YPtr() { return fY; };

    // Get W binning
    inline double* WPtr() { return fW; };
    
    // Fill array with pT values for the specified heavy-quark mass xm
    void FillPt(double* ptall, double xm);

    // Get number of contributions
    inline int GetNContr() { return fNContr; };
    
    // Get specified contribution
    inline MNRContribution* GetContr(int c) { return fContr[c]; };

    // Set cross sections in all bins to non-physical values
    void NonPhys();

    // Set cross sections in specified pT bin to non-physical values
    void NonPhys(int bpt);
    
    // Set all cross sections to zero
    void Zero();

    // Print grid (to console) for the specified heavy-quark mass
    void Print(double xm);
    
    // Find W bin that matches the specified value w
    int FindWBin(double w);

    // Transformation from original grid (gridorig) to new one (gridtrg)
    // (using cubic spline interpolation)
    static void InterpolateGrid(Grid* gridorig, Grid* gridtrg, double mq);

  // Private fields
  private:
    // Number of pT (L) bins
    int fNL;
    // Number of y bins
    int fNY;
    // Number of W bins
    int fNW;
    // pT (L) bin centres
    double* fL;
    // y bin centres
    double* fY;
    // W bin centres
    double* fW;
    // W bin boundaries
    double* fBW;
    // Cross sections in bins
    double**** fCS;
    // Number of contributions
    int fNContr;
    // Contributions
    MNRContribution** fContr;
  };
}


#endif // HVQMNR_GRID_H
