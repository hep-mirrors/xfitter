#ifndef xfitter_cpp_h
#define xfitter_cpp_h
#include "dimensions.h"
#include <string>
using namespace std;
#include <iostream>

//Fortran wrapper for lhapdferrors

extern"C" {
  void chi2_scan_();

  void mc_method_();

  //LHAPDFErrors functions
  void get_lhapdferrors_();
  void getpdfunctype_heraf_(int& lmontecarlo, int& lasymhess, int& lsymmhess, const char name[], int len);
  void set_verbosity_(int &level);
    
    
  //Error logging function
  void hf_errlog_(const int &id, const char text[], int);

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


  extern struct {
    double alpha_[NTOT_C];       // Total uncorrelated errors
    double alpha_mod_[NTOT_C];   // Total uncorrelated errors modified
    double beta_[NTOT_C][NSYSMAX_C];   // Influence of systematic errors on measurements
    double sysa_[NSYSMAX_C][NSYSMAX_C];    // Correlation matrix of systematics
    char system_[NSYSMAX_C][64];     // Names of correlated systematic errors
    int nsys_;                 // Actual number of correlated systematic sources
  } systema_;

  extern struct {
    double betaasym_[NTOT_C][2][NSYSMAX_C];      
    double omega_[NTOT_C][NSYSMAX_C];   // Quadratic term for influence of syst. errors on measurements.
    int lasymsyst_[NSYSMAX_C];        // asymmetric uncertainty
  } systasym_;

  extern struct {
    char lhapdfset_[128];
    char lhapdfvarset_[128];
    int ilhapdfset_;
    int nlhapdf_sets_;
    int nparvar_;
    int lhapdferrors_;
    int scale68_;
    int nremovepriors_;
    int lhapdfprofile_;
    int lhascaleprofile_;
  } clhapdf_;

  //aplhas
  extern struct {
    double alphas_;
  } c_alphas_;

  //boson masses
  extern struct {
    double mz_;
    double mw_;
    double mh_;
  } boson_masses_;

  extern struct ckm_matrix_cb {
    double Vud, Vus, Vub, Vcd, Vcs, Vcb, Vtd, Vts, Vtb;
  } ckm_matrix_;

  extern struct {
    int iflagfcn_;     // FCN minuit flag, set by minuit (1-init, 2-iteration, 3-finalisation)
    int nparfcn_;
    int ndfmini_;
    int ifcncount_;    // Count number of iterations
      } cfcn_;

  extern struct {
    char outdirname_[256]; // outout dir name
    char lhapdf6outdir_[256];
  } coutdirname_;

  extern struct {
    double theo_[NTOT_C];          // Theory predictions, filled for each iteration
    double theo_mod_[NTOT_C];      // Theory predictions, filled for each iteration
    double theo_fix_[NTOT_C];      // Fixed theory prediction (if given by &InTheory namelist)
    double theo_unc_[NTOT_C];      // Uncorrelated uncertainty on theory predictions 
    double theo_tot_up_[NTOT_C];	 // Total up uncertainty on theory predictions 
    double theo_tot_down_[NTOT_C]; // Total down uncertainty on theory predictions
  } c_theo_;

  extern struct {
    int npoints_;       // Actual number of data points
  } cndatapoints_;

  extern struct {
    int n_syst_meas_[NSYSMAX_C];          // Number of measurements syst. source affects
    int syst_meas_idx_[NSYSMAX_C][NTOT_C];  // data points syst. source affects	   
  } sysmeas_;

  extern struct {
    int resetcommonsyst_;
  } systematicsflags_;

  extern struct {
    double sysscalefactor_[NSYSMAX_C];
    int sysscalingtype_[NSYSMAX_C];
    int sysform_[NSYSMAX_C];
    int dooffset_;
    int chi2offsfinal_;
    int chi2offsrecalc_;
  } systscal_;

  //Penalty term for systematic sources, 1 by default.
  extern struct {
    double syspriorscale_[NSYSMAX_C];
  } csystprior_;
  
  extern struct {
    float q2val_[40];
    float starting_scale_;
    float strange_frac_;
    double chi2maxerror_;
    float hf_mass_[3];
    float charm_frac_;
    int ldebug_;
    int dobands_;
    int usegridlhapdf5_;
    int writelhapdf6_;
    int h1qcdfunc_;
    int writealphastomemberpdf_;
    int itheory_;
    int i_fit_order_;
    int iparam_;
    int hfscheme_;
    int lrand_;
    int statype_;
    int systype_;
    float outxrange_[2];
    int outnx_;
    int ilenpdf_;
    int nchebglu_;
    float chebxmin_;
    int nchebsea_;
    float wmnlen_;
    float wmxlen_;
    int ichebtypeglu_;
    int ichebtypesea_;
    int ifsttype_;
    int iseedmc_;
    int ioffsetchebsea_;
    int ewfit_;
    int npolyval_;
    int lead_;
    int izpopoly_;
    int ipolysqr_;
    int useapplg_;
    int napplgrids_;
    int lfitdy_;
    int lfastapplgrid_;
    int lranddata_;
    int idh_mod_;
    int ipdfset_;
    int icheck_qcdnum_;
    int useprevfit_;
    int corrsystbyoffset_;
    int corsysindex_;
    char statscale_[32];
    char uncorsysscale_[32];
    char corsysscale_[32];
    char uncorchi2type_[32];
    char corchi2type_[32];
    char hf_scheme_[32];
    int asymerrorsiterations_;
    int luseapplgridckm_;
    int pdfrotate_;
    int ExtraPdfs_;
    int WriteLHAPDF5_;
  } steering_;

  extern struct {
    double e_stat_poisson_[NTOT_C];
    double e_stat_const_[NTOT_C];
    double e_uncor_poisson_[NTOT_C];
    double e_uncor_const_[NTOT_C];
    double e_uncor_mult_[NTOT_C];
    double e_uncor_logNorm_[NTOT_C];
  } cuncerrors_;

  extern struct {
    char label_[64];
    double central_;
    double values_[NCHI2POINTS_C];
    int dataid_[NSET_C];
    char term_[NSET_C][16][8];
    char theorysources_[NSET_C][16][NCHI2POINTS_C][1000];
    int scan_;
    int pdferrors_;
    int pdfprofile_;
    int scaleprofile_;
    char chi2lhapdfref_[128];
    char chi2lhapdfset_[128];
    char chi2lhapdfvarset_[128];
    int chi2nparvar_;
    int chi2parpoint_;

  } chi2scan_;

  extern struct {
    int isysttype_[NSYSMAX_C];
  } csysttype_;


  extern struct {
    double men_;
    double mel_;
    double mmn_;
    double mmo_;
    double mtn_;
    double mta_;
    double mup_;
    double mdn_;
    double mch_;
    double mst_;
    double mtp_;
    double mbt_;
  } fermion_masses_;
}


#endif
