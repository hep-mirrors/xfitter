#ifndef HVQMNR_FRAG_H
#define HVQMNR_FRAG_H

#include <vector>
// Declare classes which appear in method signatures
class TH2D;
class TF1;


namespace HVQMNR
{
  //
  // Author: Oleksandr Zenaiev (oleksandr.zenaiev@desy.de)
  //
  // Class for non-perturbative heavy quark to hadron fragmentation.
  // Kartvelishvili, Peterson and BCFY fragmentation functions are provided.
  // Specific ground state mesons can be treated separately 
  // (e.g. for contribution to D^0 from D*+ decays).
  class Frag {
  // Public methods
  public:
    // Constructor
    Frag();

    // Destructor
    ~Frag();
    
    // Set number of z (rescaling variable) bins
    void SetNz(int nz);
    
    // Add final state with fragmentation function ffrag and mass M.
    // If M < 0, heavy-quark mass is used instead.
    void AddOut(TF1* ffrag, double M);
    
    // Get fragmentation function for specified final state
    TF1* GetFF(int f) { return fFFrag[f]; };
    
    // Calculate cross sections for provided grid and heavy-quark mass xm and fill histograms hcs
    void CalcCS(Grid* grid, double xm, TH2D** hcs);
    
    // Set debug flag
    void SetDebug(bool debug) { bDebug=debug; };

    // Various fragmentation functions
    // (see arXiv:hep-ph/0306212 for more details about BCFY functions and rescaling for D0 and D+ mesons)
    // BCFY for vector mesons
    static double bcfy_v(double* x, double* p);
    // BCFY for vector mesons (rescaled)
    static double bcfy_v_prim(double* x, double* p);
    // BCFY for pseudoscalar mesons
    static double bcfy_p(double* x, double* p);
    // BCFY for D0 mesons
    static double bcfy_dzero(double* x, double* p);
    // BCFY for D+ mesons
    static double bcfy_dch(double* x, double* p);
    // Kartvelishvili function [Phys.Lett. B78 (1978) 615]
    static double kar(double* x, double* p);
    // Kartvelishvili function (rescaled)
    static double kar_prim(double* x, double* p);
    // Kartvelishvili for D0 mesons
    static double kar_dzero(double* x, double* p);
    // Kartvelishvili for D+ mesons
    static double kar_dch(double* x, double* p);
    // Peterson function [Phys.Rev. D27 (1983) 105]
    static double pet(double* x, double* p);
    // Kartvelishvili with two parameters ("Misha-style" parametrisation, 
    // see DESY-THESIS-2011-033, Section 2.2 for description)
    static double karw(double* x, double* p);
    // Kartvelishvili with two parameters (rescaled)
    static double karw_prim(double* x, double* p);
    // Kartvelishvili with two parameters for D0 mesons
    static double karw_dzero(double* x, double* p);
    // Kartvelishvili with two parameters for D+ mesons
    static double karw_dch(double* x, double* p);
    // Kartvelishvili with three W bins (see arXiv:1211.1182 for description)
    static double karstep(double* x, double* p);
    // // Kartvelishvili with three W bins (rescaled)
    static double karstep_prim(double* x, double* p);
    // // Kartvelishvili with three W bins for D0 mesons
    static double karstep_dzero(double* x, double* p);
    // // Kartvelishvili with three W bins for D+ mesons
    static double karstep_dch(double* x, double* p);

    // Get fragmentation function with parameter par for specified meson.
    // Depending on ff the following fragmentaion function is returned:
    //   ff = 1:  Kartvelishvili
    //   ff = 2:  BCFY
    //   ff = 3:  Peterson
    //   ff = 10: Kartvelishvili with 2D parametrisation ("Misha-style")
    //   ff = 20: Kartvelishvili with three W bins
    // If mean pointer is provided, for 1D function it is set to mean value.
    //
    // WARNING: it is known feature that the fit might not converge if fragmentation depends 
    //          on energy in parton-parton rest frame, especially if the heavy-quark mass 
    //          is released, use ff = 10 and ff = 20 with great caution!
    static TF1* GetFragFunction(int ff, const char* meson, double par, double* mean = 0);
    
  // Private methods
  private:
    // Precalculate variables for provided grid and heavy-quark mass xm
    void Precalc(Grid* grid, double xm) ;

    // Clear z binning
    void ClearZ();
    
    // Clear precalculated variables
    void ClearPrecalc();
    
  // Public fields
  public:
    // Charm and beauty hadron masses
    static const double fM_dzero;
    static const double fM_dch;
    static const double fM_dstar;
    static const double fM_ds;
    static const double fM_lambdac;
    static const double fM_bzero;
    static const double fM_bch;
    static const double fM_bs;
    
  // Private fields
  private:
    // Final states
    int fNout; // number of final states
    std::vector<TF1*> fFFrag; // fragmentation functions
    std::vector<double> fMh2; // hadron masses squared
    
    // Precalculated variables
    int fBnz;
    double* fZc;
    double* fZw;
    std::vector<double*> fWz; 
    
    // Flags
    bool bFirst; // first run
    bool bDebug; // verbose output

    // Variables for precalculation
    double fLastxm;   // last used heavy-quark mass
    int fNRecalc;     // number of recalculations
    int fCGnpt;       // number of pT bins
    int fCGny;        // number of rapidity bins
    double* fCGpt;    // pT bins [npt]
    double* fCGy;     // rapidity bins [ny]
    double** fCGptf;  // contribution in rescaled pT bins [npt][nz]
    double**** fCGyf; // contribution in rescaled pT-rapidity bins [npt][ny][nz][nf]
  };
}


#endif // HVQMNR_FRAG_H
