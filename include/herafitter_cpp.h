#ifndef herafitter_cpp_h
#define herafitter_cpp_h

using namespace std;

//Fortran wrapper for lhapdferrors

extern"C" {
  //LHAPDFErrors functions
  void get_lhapdferrors_();
  void getpdfunctype_heraf_(const char name[30], int& lmontecarlo, int& lasymhess, int& lsymmhess);

  //Error logging function
  void hf_errlog_(const int &id, const char text[], int);

  double chi2data_theory_(const int &iflag);

  //IO functions
  void writefittedpoints_();
  void store_pdfs_(const char filename[], int);
  void writetheoryfiles_(const int& nnuisance, double theo_cent[2500], const int& symmetricpdferr);
  void fopen_(const int& fnumber, const char fname[], int);
  void fclose_(const int& fnumber);

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
			     int& Ncorrelated, double* Uncor);


  extern struct {
    double alpha_[2500];       // Total uncorrelated errors
    double alpha_mod_[2500];   // Total uncorrelated errors modified
    double beta_[2500][300];   // Influence of systematic errors on measurements
    double sysa_[300][300];    // Correlation matrix of systematics
    char system_[300][64];     // Names of correlated systematic errors
    int nsys_;                 // Actual number of correlated systematic sources
  } systema_;

  extern struct {
    double betaasym_[2500][2][300];      
    double omega_[2500][300];   // Quadratic term for influence of syst. errors on measurements.
    int lasymsyst_[300];        // asymmetric uncertainty
  } systasym_;

  extern struct {
    char lhapdfset_[128];
    char lhapdfvarset_[128];
    int ilhapdfset_;
    int nlhapdf_sets_;
    int nparvar_;
    int lhapdferrors_;
    int scale68_;
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

  extern struct {
    int iflagfcn_;     // FCN minuit flag, set by minuit (1-init, 2-iteration, 3-finalisation)
    int nparfcn_;
    int ndfmini_;
    int ifcncount_;    // Count number of iterations
      } cfcn_;

  extern struct {
    char outdirname_[128]; // outout dir name
    char lhapdf6outdir_[128];
  } coutdirname_;

  extern struct {
    double theo_[2500];          // Theory predictions, filled for each iteration
    double theo_mod_[2500];      // Theory predictions, filled for each iteration
    double theo_fix_[2500];      // Fixed theory prediction (if given by &InTheory namelist)
    double theo_unc_[2500];      // Uncorrelated uncertainty on theory predictions 
    double theo_tot_up_[2500];	 // Total up uncertainty on theory predictions 
    double theo_tot_down_[2500]; // Total down uncertainty on theory predictions
  } c_theo_;

  extern struct {
    int npoints_;       // Actual number of data points
  } cndatapoints_;

  extern struct {
    int n_syst_meas_[300];          // Number of measurements syst. source affects
    int syst_meas_idx_[300][2500];  // data points syst. source affects	   
  } sysmeas_;

  extern struct {
    int resetcommonsyst_;
  } systematicsflags_;

  extern struct {
    double sysscalefactor_[300];
    int sysscalingtype_[300];
    int sysform_[300];
    int dooffset_;
    int chi2offsfinal_;
    int chi2offsrecalc_;
  } systscal_;

  /*
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
    real wmnlen_;
    real wmxlen_;
    int ichebtypeglu_;
    int ichebtypesea_;
    int ifsttype_;
    int iseedmc_;
      common/STEERING/
     $     ,IOFFSETCHEBSEA,EWFIT
     $     ,npolyval, lead
     $     ,IZPOPOLY,IPOLYSQR, useapplg, napplgrids,lfitdy
     $     ,LFastAPPLGRID,Lranddata,iDH_MOD
     $     ,IPDFSET, ICHECK_QCDNUM
     $     ,UsePrevFit
     $     ,CorrSystByOffset,CorSysIndex
     $     ,StatScale, UncorSysScale, CorSysScale,UncorChi2Type,
     $     CorChi2Type, hf_scheme, AsymErrorsIterations,
     $     LUseAPPLgridCKM
  } steering_;
  */

  extern struct {
    double e_stat_poisson_[2500];
    double e_stat_const_[2500];
    double e_uncor_poisson_[2500];
    double e_uncor_const_[2500];
    double e_uncor_mult_[2500];
    double e_uncor_logNorm_[2500];
  } cuncerrors_;
}

#endif
