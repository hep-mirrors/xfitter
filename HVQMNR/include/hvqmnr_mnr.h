#ifndef HVQMNR_MNR_H
#define HVQMNR_MNR_H


namespace HVQMNR
{
  //
  // Author: Oleksandr Zenaiev (oleksandr.zenaiev@desy.de)
  //
  // Class for MNR calculation: heavy-quark production in hadron collisions at NLO
  // (using original FORTRAN routines from hvqcrsx.f)
  // [M. Mangano, P. Nason and G. Ridolfi, Nucl. Phys. B 373 (1992) 295.]
  class MNR {

  // Public members
  public:
    // Constructor
    MNR();
    
    // Destructor
    ~MNR();
    
    // Set perturbative scale coefficients
    //
    // Scales are parametrised as:
    // mu_f^2 = mf_a * pT^2 + mf_b * xm^2 + mf_c
    // mu_r^2 = mr_a * pT^2 + mr_b * xm^2 + mr_c
    // where mu_f, mu_r are factorisation and renormalisation, respectively,
    // pT is transverse momentum and xm is the heavy-quark mass.
    void SetScaleCoef(double mf_a, double mf_b, double mf_c, double mr_a, double mr_b, double mr_c);
    
    // Set debug flag
    void SetDebug(int debug) { bDebug = debug; };

    // Calculate constants
    void CalcConstants();
    
    // Calculate binning
    void CalcBinning();

    // Calculate cross sections for provided grid and heavy-quark mass xm
    void CalcXS(Grid* grid, double xm);

        
  // Private members
  private:
    // Get factorisation scale
    double GetMf2(double xm2, double pt2);
    // Get renormalisation scale
    double GetMr2(double xm2, double pt2);

    // Precalculate PDFs at the specified factorisation scale
    void PrecalculatePDF(double mf2);
    // Get precalculated PDFs in form of needed "structure functions" (gg, qq and qg)
    int GetSF(double& pdf_gg, double& pdf_qq, double& pdf_qq_a, double& pdf_qg, double& pdf_qg_a, double& pdf_qg_r, double& pdf_qg_a_r, double adoptx1, double adoptx2, double mf2 = -1.0);

    // Get PDFs (recalculate, e.g. call to xFitter routine)
    void GetPDF(double mf2, double x, double pdf[13]);
    // Get strong coupling
    double GetAs(double mr2);

    // Perform precalculation
    void Precalc(Grid* grid);

  // Public fields
  public:
    // Centre-of-mass energy squared
    double fC_sh;      

    // Number of light flavours
    int fC_nl;

    // Number of y bins in centre-of-mass system
    int fBn_x3;

    // Number of t3 bins (3 body variable)
    int fBn_x4;

    // Number of points in x used for fast structure functions evaluation
    int fSF_nb;

    // Contrbution flags
    bool bFS_Q; // particle final state
    bool bFS_A; // antiparticle final state
    
  // Private fields
  private:
    // Constants
    const static double fC_pi;
    const static double fC_2pi;
    const static double fC_hc2;
    const static double fC_vca;
    const static double fC_vtf;
    double fC_b0;
    // Centre-of-mass energy
    double fC_sqrt_sh;    
    // Normalisation factor
    double fC_xnorm;
    // Perturbative scale coefficients
    double fMf_A;
    double fMf_B;
    double fMf_C;
    double fMr_A;
    double fMr_B;
    double fMr_C;

    // Variables for fast structure functions evaluation
    const static int fSF_npart;
    const static double fSF_min_x;
    const static double fSF_max_x;
    const static double fSF_min_mf2;
    const static double fSF_max_mf2;
    double fSF_log10_min_x;
    double fSF_log10_max_x;
    double fSF_min_adoptx;
    double fSF_max_adoptx;
    double fSF_step_log10_x;
    double* fSF_pdf;

    // Rapidity bins in centre-of-mass system
    double* fBc_x3;
    double* fBw_x3;

    // t3 bins (3 body variable)
    double* fBc_x4;
    double* fBw_x4;
    
    // Precalcuated grid variables
    int fNRecalc;
    // Pointer to all allocated memory
    double* fC_mem;
    // LO gg
    double* fCh0_hqh0gg;
    // LO qq
    double* fCh0_hqh0qa;
    // NLO gg
    double* fCh3_hqhpgg;
    double* fCh3_hqbpgg;
    double* fCh3_hqhlgg;
    double* fCh2_hqhdgg;
    double* fCh2_hqbdgg;
    double* fCh2_hqh0gg;
    double* fCh3c_hqhpgg;
    double* fCh3c_hqbpgg;
    double* fCh3c_hqhlgg;
    // NLO qq
    double* fCh3_hqhpqa;
    double* fCh3_hqbpqa;
    double* fCh3_hqhlqa;
    double* fCh3_a_ashpqa;
    double* fCh2_hqhdqa;
    double* fCh2_hqbdqa;
    double* fCh2_hqh0qa;
    double* fCh2_a_ashdqa;
    double* fCh3c_hqhpqa;
    double* fCh3c_hqbpqa;
    double* fCh3c_hqhlqa;
    double* fCh3c_a_ashpqa;
    // NLO qg
    double* fCh3_hqhpqg;
    double* fCh3_hqbpqg;
    double* fCh3_hqhlqg;
    double* fCh3_r_hqhpqg;
    double* fCh3_r_hqbpqg;
    double* fCh3_r_hqhlqg;
    double* fCh3_a_ashpqg;
    double* fCh3_a_r_ashpqg;
    // Normalisation coefficients
    double* fC_N;
    double* fC_NN;
    // Kinematics
    double* fCk_t32;//[sBn_x4]
    double* fCk_t33;//[sBn_x4]
    double* fCk_tx;//[sBn_x4]
    double* fCk_lntx;//[sBn_x4]
    double* fCk_lntx_o_tx;//[sBn_x4]
    double* fCk_pxtcor;//[sBn_x4]
    double fCk_sum_o_tx;
    double fCk_sum_lntx_o_tx;
    double* fCk_t1t;//[sBn_x3]
    double* fCk_chyprim2;//[sBn_x3]

    // Flags
    bool bFirst; // first run
    bool bDebug; // verbose output
  };
}


#endif // HVQMNR_MNR_H
