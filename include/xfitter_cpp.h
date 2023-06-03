#ifndef xfitter_cpp_h
#define xfitter_cpp_h

#include "dimensions.h"
#include <string>
using namespace std;
#include <iostream>
#include "xfitter_cpp_base.h"

//Fortran wrapper for lhapdferrors

extern"C" {
  void chi2_scan_();

  void mc_method_();

  //LHAPDFErrors functions
  void get_lhapdferrors_();
  void getpdfunctype_heraf_(int& lmontecarlo, int& lasymhess, int& lsymmhess, const char name[], int len);
  void set_verbosity_(int &level);
    
    
  double chi2data_theory_(const int &iflag);

  //IO functions
  void writefittedpoints_();
  void store_pdfs_(const char filename[], int);
  void writetheoryfiles_(const int& nnuisance, double theo_cent[NTOT_C], const int& symmetricpdferr);
  void fopen_(const int& fnumber, const char fname[], int);
  void fclose_(const int& fnumber);

  //banner
  void hfbanner_();

  //lhapdf6 functions
  void fill_c_common_();
  void print_lhapdf6_();
  void save_data_lhapdf6_(int& iset);

  //applgrid pdf and alphas functions: externally defined alpha_s and pdf routines for fortran callable convolution wrapper
  void appl_fnpdf_(const double& x, const double& Q, double* f);
  double appl_fnalphas_(const double& Q);

  //Covariance matrix to nuisance parameter conversion
  void getnuisancefromcovar_(const int& NDimCovar, const int& NDimSyst, const int& NCovar,
			     double* Covar, double *ANuisance, const double& Tolerance, 
			     int& Ncorrelated, double* Uncor, const int& LSepDiag);

  //chi2 evaluation
  void getnewchisquare_(int &flag_in, int &n0_in, double &fchi2_in, double *rsys_in, double *ersys_in, double *pchi2_in, double &fcorchi2_in);

  extern struct {
    double alpha[NTOT_C];       // Total uncorrelated errors
    double alpha_mod[NTOT_C];   // Total uncorrelated errors modified
    double beta[NTOT_C][NSYSMAX_C];   // Influence of systematic errors on measurements
    double sysa[NSYSMAX_C][NSYSMAX_C];    // Correlation matrix of systematics
    char system[NSYSMAX_C][64];     // Names of correlated systematic errors
    int nsys;                 // Actual number of correlated systematic sources
  } systema_;

  extern struct {
    double betaasym[NTOT_C][2][NSYSMAX_C];      
    double omega[NTOT_C][NSYSMAX_C];   // Quadratic term for influence of syst. errors on measurements.
    int lasymsyst[NSYSMAX_C];        // asymmetric uncertainty
  } systasym_;

  extern struct {
    double scgamma_[NTOT_C][NSYSMAX_C]; //scaled gamma
    double scomega_[NTOT_C][NSYSMAX_C]; //scaled omega
    double sysshift_[NSYSMAX_C];        //systematic shift
    double scerrors_[NTOT_C];           //scaled uncorrelated errors
  } systexport_;
  
  extern struct {
    double daten_[NTOT_C];   //!> Data values
    double e_unc_[NTOT_C];              //!> Uncorelated uncertainty 
    double e_tot_[NTOT_C];              //!> Total uncertainty
    double e_sta_[NTOT_C];              //!> Uncorelated uncertainty
    double e_sta_const_[NTOT_C];              //!> Uncorelated uncertainty 
    double e_unc_const_[NTOT_C];        //!> Uncorelated, unscaled error
    int jset_[NTOT_C];        //!> Uncorelated, unscaled error
    int indextheorybin_[NTOT_C];
    int jplot_[NTOT_C];
  } indata2_;
  
  extern struct {
    char lhapdfset[128];
    char lhapdfvarset[128];
    int ilhapdfset;
    int nlhapdf_sets;
    int nparvar;
    int lhapdferrors;
    int scale68;
    int nremovepriors;
    int lhapdfprofile;
    int lhascaleprofile;
  } clhapdf_;

  //alphas
  extern struct {
    double alphas;
  } c_alphas_;

  //boson masses
  extern struct {
    double Mz;
    double Mw;
    double Mh;
  } boson_masses_;

  // Widths:
  extern struct {
    double Wz, Ww, Wh, Wtp;
  } widths_;

  // CKM
  extern struct ckm_matrix_cb {
    double Vud, Vus, Vub, Vcd, Vcs, Vcb, Vtd, Vts, Vtb;
  } ckm_matrix_;

  // EW coupligs:
  extern struct {
    double Alphaem;
    double sin2thW;
    double cos2thW;
  } ew_couplings_;

    // Constants:
  extern struct {
    double Gf, ConvFac;
  } constants_;


  extern struct {
    int iflagfcn;     // FCN minuit flag, set by minuit (1-init, 2-iteration, 3-finalisation)
    int nparfcn;
    int ndfmini;
    int ifcncount;    // Count number of iterations
      } cfcn_;

  extern struct {
    char outdirname[256]; // output dir name
  } coutdirname_;

  extern struct {
    double theo[NTOT_C];          // Theory predictions, filled for each iteration
    double theo_mod[NTOT_C];      // Theory predictions, filled for each iteration
    double theo_fix[NTOT_C];      // Fixed theory prediction (if given by &InTheory namelist)
    double theo_unc[NTOT_C];      // Uncorrelated uncertainty on theory predictions 
    double theo_tot_up[NTOT_C];	 // Total up uncertainty on theory predictions 
    double theo_tot_down[NTOT_C]; // Total down uncertainty on theory predictions
  } c_theo_;

  extern struct {
    int npoints;       // Actual number of data points
  } cndatapoints_;

  extern struct {
    int n_syst_meas[NSYSMAX_C];          // Number of measurements syst. source affects
    int syst_meas_idx[NSYSMAX_C][NTOT_C];  // data points syst. source affects	   
  } sysmeas_;

  extern struct {
    int resetcommonsyst;
  } systematicsflags_;

  extern struct {
    double sysscalefactor[NSYSMAX_C];
    int sysscalingtype[NSYSMAX_C];
    int sysform[NSYSMAX_C];
    int dooffset;
    int chi2offsfinal;
    int chi2offsrecalc;
  } systscal_;

  //Penalty term for systematic sources, 1 by default.
  extern struct {
    double syspriorscale[NSYSMAX_C];
  } csystprior_;
  
  extern struct {
    float q2val[40];
    double chi2maxerror;
    int ldebug;
    int i_fit_order;
    int lrand;
    int statype;
    int systype;
    float outxrange[2];
    int outnx;
    int iseedmc;
    int useapplg;
    int napplgrids;
    int lranddata;
    int idh_mod;
    int ipdfset;
    int vIPDFSET;
    int cIDataSet;
    int icheck_qcdnum;
    int useprevfit;
    int corrsystbyoffset;
    int corsysindex;
    char statscale[32];
    char uncorsysscale[32];
    char corsysscale[32];
    char uncorchi2type[32];
    char corchi2type[32];
    int asymerrorsiterations;
    int pdfrotate;
    int ExtraPdfs;
    int WriteLHAPDF5;
    int steering_check;   // Keep this always last
  } steering_;

  extern struct {
    double e_stat_poisson[NTOT_C];
    double e_stat_const[NTOT_C];
    double e_uncor_poisson[NTOT_C];
    double e_uncor_const[NTOT_C];
    double e_uncor_mult[NTOT_C];
    double e_uncor_logNorm[NTOT_C];
  } cuncerrors_;

  extern struct {
    char label[64];
    double central;
    double values[NCHI2POINTS_C];
    int dataid[NSET_C];
    char term[NSET_C][NTERMSMAX_C][8];
    char theorysources[NSET_C][NTERMSMAX_C][NCHI2POINTS_C][1000];
    int scan;
    int pdferrors;
    int pdfprofile;
    int scaleprofile;
    char chi2lhapdfref[128];
    char chi2lhapdfset[128];
    char chi2lhapdfvarset[128];
    int chi2nparvar;
    int chi2parpoint;

  } chi2scan_;

  extern struct {
    int isysttype[NSYSMAX_C];
  } csysttype_;


  extern struct {
    double men;
    double mel;
    double mmn;
    double mmo;
    double mtn;
    double mta;
    double mup;
    double mdn;
    double mch;
    double mst;
    double mtp;
    double mbt;
  } fermion_masses_;


  extern struct {
    double residuals_[NTOT_C];
  } c_resid_;

  extern struct {
    double chi2_poi_tot;
    double chi2_poi[NSET_C];
    double chi2_poi_data[NTOT_C];
  } cdatapoi_;

  extern struct {
    int chi2poissoncorr;
    int chi2firstiterationrescale;
    int chi2extrasystrescale;
  } chi2options_;
}

#endif
