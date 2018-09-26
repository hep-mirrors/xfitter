// emacs: this is -*- c++ -*-

//  appl_grid.h       

//  grid class header - all the functions needed to create and 
//  fill the grid from an NLO calculation program
//  
//  Copyright (C) 2007 Mark Sutton (sutt@hep.ucl.ac.uk)    

// $Id: appl_grid.h, v1.00 2007/10/16 17:01:39 sutt

// Fixme: this needs to be tidied up. eg there are too many different, 
//        and too many version of, accessors for x/y, Q2/tau etc there 
//        should be only  one set, for x and Q2 *or* y and tau, but 
//        not both. In fact they should be for x and Q2, since y and tau 
//        should be purely an internal grid issue of no concern for the 
//        user.

#ifndef __APPL_GRID_H
#define __APPL_GRID_H

#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <exception>

#include "TH1D.h"


double _fy(double x);
double _fx(double y);
double _fun(double y);


#include "correction.h"

namespace appl { 


/// forward declarations - full definitions included 
/// from appl_grid.cxx 
class igrid;
class appl_pdf;


const int MAXGRIDS = 5;


/// externally visible grid class
class grid {

public:

  // grid error exception
  class exception : public std::exception { 
  public:
    exception(const std::string& s) { std::cerr << what() << " " << s << std::endl; }; 
    //exception(std::ostream& s)      { std::cerr << what() << " " << s << std::endl; }; 
    exception(std::ostream& s)      { std::stringstream ss; ss << s.rdbuf(); std::cerr << what() << " " << ss.str() << std::endl; }; 
    virtual const char* what() const throw() { return "appl::grid::exception"; }
  };

  typedef enum { STANDARD=0, AMCATNLO=1, SHERPA=2, LAST_TYPE=3 } CALCULATION; 

public:

  grid(int NQ2=50,  double Q2min=10000.0, double Q2max=25000000.0,  int Q2order=5,  
       int Nx=50,   double xmin=1e-5,     double xmax=0.9,          int xorder=5,
       int Nobs=20, double obsmin=100.0,  double obsmax=7000.0, 
       std::string genpdf="mcfm_pdf", 
       int leading_order=0, int nloops=1, 
       std::string transform="f2");

  grid( int Nobs, const double* obsbins,
	int NQ2=50,  double Q2min=10000.0, double Q2max=25000000.0, int Q2order=5,
        int Nx=50,   double xmin=1e-5,     double xmax=0.9,         int xorder=5, 
	std::string genpdf="mcfm_pdf",
	int leading_order=0, int nloops=1, 
	std::string transform="f2" );

  grid( const std::vector<double>& obs,
	int NQ2=50,  double Q2min=10000.0, double Q2max=25000000.0,   int Q2order=5, 
        int Nx=50,   double xmin=1e-5,     double xmax=0.9,           int xorder=5, 
	std::string genpdf="mcfm_pdf", 
	int leading_order=0, int nloops=1, 
	std::string transform="f2" );

  // build a grid but don't build the internal igrids - these can be added later
  grid( const std::vector<double>& obs,
	std::string genpdf="nlojet_pdf", 
	int leading_order=0, int nloops=1, 
	std::string transform="f2" );

  // copy constructor
  grid(const grid& g);

  // read from a file
  grid(const std::string& filename="./grid.root", const std::string& dirname="grid");

  // add an igrid for a given bin and a given order 
  void add_igrid(int bin, int order, igrid* g);

  virtual ~grid();
  
  // update grid with one set of event weights
  void fill(const double x1, const double x2, const double Q2, 
	    const double obs, 
	    const double* weight, const int iorder);
  
  
  void fill_phasespace(const double x1, const double x2, const double Q2, 
		       const double obs, 
		       const double* weight, const int iorder);
  
  
  void fill_grid(const double x1, const double x2, const double Q2, 
		 const double obs, 
		 const double* weight, const int iorder)  {
    if (isOptimised())   fill(x1, x2, Q2, obs, weight, iorder);
    else                 fill_phasespace(x1, x2, Q2, obs, weight, iorder);
  }


  void fill_index(const int ix1, const int ix2, const int iQ2, 
		  const int iobs, 
		  const double* weight, const int iorder);


  // trim/untrim the grid to reduce memory footprint
  void trim();
  void untrim();
 
  // formatted output 
  std::ostream& print(std::ostream& s=std::cout) const;

  // don't do anything anymore
  // void setuppdf(void (*pdf)(const double& , const double&, double* ) );

  // get the interpolated pdf's
  //  void pdfinterp(double x1, double Q2, double* f);


  // perform the convolution to a specified number of loops
  // nloops=-1 gives the nlo part only
  std::vector<double>  vconvolute(void   (*pdf)(const double& , const double&, double* ), 
				  double (*alphas)(const double& ), 
				  int     nloops, 
				  double  rscale_factor=1,
				  double  fscale_factor=1,
				  double  Escale=1 );

  std::vector<double>  vconvolute(void   (*pdf1)(const double& , const double&, double* ), 
				  void   (*pdf2)(const double& , const double&, double* ), 
				  double (*alphas)(const double& ), 
				  int     nloops, 
				  double  rscale_factor=1,
				  double  fscale_factor=1,
				  double  Escale=1 );


  // perform the convolution to a specified number of loops
  // nloops=-1 gives the nlo part only
  std::vector<double>  vconvolute(double Escale,
				  void   (*pdf)(const double& , const double&, double* ), 
				  double (*alphas)(const double& ), 
				  int     nloops, 
				  double  rscale_factor=1,
				  double  fscale_factor=1  ) { 
    return vconvolute(pdf, alphas, nloops, rscale_factor, fscale_factor, Escale); 
  }


  // perform the convolution to the max number of loops in grid
  std::vector<double> vconvolute(void   (*pdf)(const double& , const double&, double* ), 
				 double (*alphas)(const double& ) )   { 
    return vconvolute( pdf, alphas, m_order-1 ); 
  } 


  // perform the convolution to the max number of loops in grid
  std::vector<double> vconvolute(double Escale, 
				 void   (*pdf)(const double& , const double&, double* ), 
				 double (*alphas)(const double& ) )   { 
    return vconvolute( Escale, pdf, alphas, m_order-1 ); 
  } 


  // perform the convolution to a specified number of loops 
  // for a single sub process, nloops=-1 gives the nlo part only
  std::vector<double> vconvolute_subproc(int subproc, 
					 void   (*pdf)(const double& , const double&, double* ), 
					 double (*alphas)(const double& ), 
					 int     nloops, 
					 double  rscale_factor=1, double Escale=1 ); 


  // perform the convolution to a specified number of loops 
  // for a single sub process, nloops=-1 gives the nlo part only
  std::vector<double> vconvolute_subproc(int subproc, 
					 double Escale, 
					 void   (*pdf)(const double& , const double&, double* ), 
					 double (*alphas)(const double& ), 
					 int     nloops, 
					 double  rscale_factor=1 ) { 
    return vconvolute_subproc(subproc, pdf, alphas, nloops, rscale_factor, Escale); 
  } 


  // perform the convolution to the max number of loops in grid 
  // for a single sub process
  std::vector<double> vconvolute_subproc(int subproc, 
					 void   (*pdf)(const double& , const double&, double* ), 
					 double (*alphas)(const double& ) )   { 
    return vconvolute_subproc( subproc, pdf, alphas, m_order-1 ); 
  } 

  // perform the convolution to the max number of loops in grid 
  // for a single sub process
  std::vector<double> vconvolute_subproc(int subproc, 
					 double Escale,
					 void   (*pdf)(const double& , const double&, double* ), 
					 double (*alphas)(const double& ) )   { 
    return vconvolute_subproc( subproc, Escale, pdf, alphas, m_order-1 ); 
  } 
  

  double vconvolute_bin( int bin, 
			 void (*pdf)(const double& , const double&, double* ), 
			 double (*alphas)(const double&) ); 
				       

  // perform the convolution to a specified number of loops
  // nloops=-1 gives the nlo part only
  TH1D* convolute(void   (*pdf)(const double& , const double&, double* ), 
		  double (*alphas)(const double& ), 
		  int     nloops, 
		  double  rscale_factor=1,
		  double  fscale_factor=1,
		  double  Escale=1 );

  // perform the convolution to a specified number of loops
  // nloops=-1 gives the nlo part only
  TH1D* convolute(void   (*pdf1)(const double& , const double&, double* ), 
		  void   (*pdf2)(const double& , const double&, double* ), 
		  double (*alphas)(const double& ), 
		  int     nloops, 
		  double  rscale_factor=1,
		  double  fscale_factor=1,
		  double  Escale=1 );


  TH1D* convolute(double Escale,
		  void   (*pdf)(const double& , const double&, double* ), 
		  double (*alphas)(const double& ), 
		  int     nloops, 
		  double  rscale_factor=1,
		  double  fscale_factor=1 ) { 
    return convolute(pdf, alphas, nloops, rscale_factor, fscale_factor, Escale); 
  }


  // perform the convolution to the max number of loops in grid
  TH1D* convolute(void   (*pdf)(const double& , const double&, double* ), 
		  double (*alphas)(const double& ) )   { 
    return convolute( pdf, alphas, m_order-1 ); 
  } 

  // perform the convolution to the max number of loops in grid
  TH1D* convolute(double Escale,
		  void   (*pdf)(const double& , const double&, double* ), 
		  double (*alphas)(const double& ) )   { 
    return convolute( Escale, pdf, alphas, m_order-1 ); 
  } 


  // perform the convolution to a specified number of loops 
  // for a single sub process, nloops=-1 gives the nlo part only
  TH1D* convolute_subproc(int subproc, 
			  void   (*pdf)(const double& , const double&, double* ), 
			  double (*alphas)(const double& ), 
			  int     nloops, 
			  double  rscale_factor=1, double Escale=1 );

  TH1D* convolute_subproc(int subproc, 
			  double Escale,
			  void   (*pdf)(const double& , const double&, double* ), 
			  double (*alphas)(const double& ), 
			  int     nloops, 
			  double  rscale_factor=1 ) { 
    return convolute_subproc( subproc, pdf, alphas, nloops, rscale_factor, Escale);
  }

  // perform the convolution to the max number of loops in grid 
  // for a single sub process
  TH1D* convolute_subproc(int subproc, 
			  void   (*pdf)(const double& , const double&, double* ), 
			  double (*alphas)(const double& ) )   { 
    return convolute_subproc( subproc, pdf, alphas, m_order-1 ); 
  } 

  TH1D* convolute_subproc(int subproc, 
			  double Escale,
			  void   (*pdf)(const double& , const double&, double* ), 
			  double (*alphas)(const double& ) )   { 
    return convolute_subproc( subproc, Escale, pdf, alphas, m_order-1 ); 
  } 


  // optimise the bin limits
  void optimise(bool force=false);
  void optimise(int NQ2, int Nx);
  void optimise(int NQ2, int Nx1, int Nx2);

  // redefine the limits by hand
  void redefine(int iobs, int iorder,
		int NQ2, double Q2min, double Q2max, 
		int Nx,  double  xmin, double  xmax);

  bool setNormalised(bool t=true) { return m_normalised=t; } 
  bool getNormalised() const      { return m_normalised; } 


  // set the filling to be symmetric and test status
  bool symmetrise(bool t=true) { return m_symmetrise=t; } 
  bool isSymmetric()     const { return m_symmetrise; } 

  bool reweight(bool t=false); 

  // access to internal grids if need be
  const igrid* weightgrid(int iorder, int iobs) const { return m_grids[iorder][iobs]; }
  
  // save grid to specified file
  void Write(const std::string& filename, const std::string& dirname="grid", const std::string& pdfname="" );

  // accessors for the observable after possible bin combination
  int    Nobs()               const { return m_obs_bins_combined->GetNbinsX(); }
  double obs(int iobs)        const { return m_obs_bins_combined->GetBinCenter(iobs+1); } 
  int    obsbin(double obs)   const { return m_obs_bins_combined->FindBin(obs)-1; } 
  double obslow(int iobs)     const { return m_obs_bins_combined->GetBinLowEdge(iobs+1); }
  double obsmin()             const { return obslow(0); } 
  double obsmax()             const { return obslow(Nobs()); } 
  double deltaobs(int iobs)   const { return m_obs_bins_combined->GetBinWidth(iobs+1); }

  const TH1D* getReference() const { return m_obs_bins_combined; } 
  TH1D*       getReference()       { return m_obs_bins_combined; } 


  //  TH1D*  getXReference() {
  //    combineReference(); 
  //    return m_obs_bins_combined; 
  //  }
 
 
  // accessors for the observable befor any bin combination
  int    Nobs_internal()               const { return m_obs_bins->GetNbinsX(); }
  double obs_internal(int iobs)        const { return m_obs_bins->GetBinCenter(iobs+1); } 
  int    obsbin_internal(double obs)   const { return m_obs_bins->FindBin(obs)-1; } 
  double obslow_internal(int iobs)     const { return m_obs_bins->GetBinLowEdge(iobs+1); }
  double deltaobs_internal(int iobs)   const { return m_obs_bins->GetBinWidth(iobs+1); }
  double obsmin_internal()             const { return obslow_internal(0); } 
  double obsmax_internal()             const { return obslow_internal(Nobs_internal()); } 

  const TH1D* getReference_internal() const { return m_obs_bins; } 
  TH1D*       getReference_internal()       { return m_obs_bins; } 

 
 

  // number of subprocesses 
  int subProcesses(int i) const;

  // general status accessors
  double& run() { return m_run; }
 
  // accessors for the status information
  bool isOptimised() const { return m_optimised; }
  bool isTrimmed()   const { return m_trimmed; }

  // lowest order of process
  int  leadingOrder() const { return m_leading_order; } 

  /// maximum number of orders ( lo=1, nlo=2, nnlo=3 )  
  /// but aMC@NLO uses 4 grids for the NLO, so m_order
  /// will be 4, but really it is still only available 
  /// 1 loop, so take account of this 
  int  nloops() const { 
    if ( m_type!=AMCATNLO ) return m_order-1;
    else if ( m_order>0 ) return 1;
    else                  return 0;
  } 

  // find out which transform and which pdf combination are being used
  std::string getTransform() const { return m_transform; }

  static double transformvar();
  static double transformvar(double v);

  std::string getGenpdf()    const { return m_genpdfname; }

  std::string version()      const { return m_version; } 
  std::string appl_version() const;

  double getCMSScale()          const { return m_cmsScale; }
  void   setCMSScale(double cmsScale) { m_cmsScale=cmsScale; }

  double getDynamicScale()          const     { return m_dynamicScale; }
  void   setDynamicScale(double dynamicScale) { m_dynamicScale=dynamicScale; }


  // set optimise flag on all sub grids
  bool setOptimised(bool t=true) { 
    return m_optimised=t;
    //    for ( int iorder=0 ; iorder<2 ; iorder++ ) { 
    //      for ( int iobs=0 ; iobs<Nobs() ; iobs++ ) m_grids[iorder][iobs]->setOptimised(t); 
    //    }
  }

  // find the number of words used for storage
  int size() const; 

  // get the cross sections
  double& crossSection()      { return m_total; } 
  double& crossSectionError() { return m_totalerror; } 
  
  //  double Lambda() const { return m_Lambda2; }

  // very lovely algebraic operators
  grid& operator=(const grid& g); 
  grid& operator*=(const double& d); 
  grid& operator+=(const grid& g);

  /// test if grids have the same limits etc
  bool operator==(const grid& g) const;   

  // shouldn't have these, the grid is too large a structure 
  // to be passed in a return
  // grid operator*(const double& d) const { return grid(*this)*=d; }
  // grid operator+(const grid& g)   const { return grid(*this)+=g; }
  
  void setDocumentation(const std::string& s);
  void addDocumentation(const std::string& s);

  std::string  getDocumentation() const { return m_documentation; }
  std::string& getDocumentation()       { return m_documentation; }


  /// set the range of the observable bins, with an optional
  /// scaling of the observable valuesfor channging units
  void setBinRange(int ilower, int iupper, double xScaleFactor=1);
  void setRange(double lower, double upper, double xScaleFactor=1);


  /// add a correction as a std::vector
  void addCorrection( std::vector<double>& v, const std::string& label="", bool combine=false );


  /// add a correction by histogram
  void addCorrection(TH1D* h, const std::string& label="", double scale=1, bool combine=false );

  
  /// access the corrections
  //  const std::vector<std::vector<double> >& corrections() const { 
  const std::vector<correction>& corrections() const { 
    return m_corrections;
  }

  /// get the correction labels
  const std::vector<std::string >& correctionLabels() const { 
    return m_correctionLabels;
  }

  /// will the corrections be applied? 
  bool getApplyCorrections() const { return m_applyCorrections; } 

  bool setApplyCorrections(bool b) { 
    std::cout << "appl::grid bin-by-bin corrections will " 
	      << ( b ? "" : "not " ) << "be applied" << std::endl;
    return m_applyCorrections=b; 
  } 

  /// apply corrections to a std::vector
  void applyCorrections(std::vector<double>& v, std::vector<bool>& applied);


  /// will a specific correction be applied? 
  bool getApplyCorrection(unsigned i) const {
    if      ( m_applyCorrections )         return true; 
    else if ( i<m_applyCorrection.size() ) return m_applyCorrection.at(i);
    return false; 
  }
 
  bool setApplyCorrection(unsigned i, bool b) { 
    if ( i>=m_corrections.size() ) return false; 
    std::cout << "appl::grid bin-by-bin correction will " 
	      << ( b ? "" : "not " ) << "be applied for correction " << i;
    if ( m_correctionLabels[i]!="" ) std::cout << " ("  << m_correctionLabels[i] << ")";
    std::cout << std::endl;
    return m_applyCorrection[i]=b;
  } 
  
  /// apply corrections to a std::vector
  bool applyCorrection(unsigned i, std::vector<double>& v);
  

  /// set the ckm matrix values if need be
  /// takes a 3x3 matrix with the format { { Vud, Vus, Vub }, { Vcd, Vcs, Vcb }, { Vtd, Vts, Vtb } }  
  void setckm( const std::vector<std::vector<double> >& ckm );

  /// takes a flat 9 element vector (or c array) with the format { Vud, Vus, Vub, Vcd, Vcs, Vcb, Vtd, Vts, Vtb }  
  void setckm( const std::vector<double>& ckm );
  void setckm( const double* ckm );


  /// set the squared ckm matrix values if need be
  /// the squared terms for eihter W+ or W- production - you probably should use setckm()
  void setckm2( const std::vector<std::vector<double> >& ckm2 );

  /// set the ckm matrix and squared ckm matrix values if need be
  const std::vector<std::vector<double> >& getckm()  const;
  const std::vector<std::vector<double> >& getckm2() const;


  /// flag custom convolution routines

  void sherpa()   { m_type = SHERPA;   std::cout << "appl::grid::sherpa()   using SHERPA convolution" << std::endl; }
  void amcatnlo() { m_type = AMCATNLO; std::cout << "appl::grid::amcatnlo() using aMC@NLO convolution" << std::endl; }
  void standard() { m_type = STANDARD; std::cout << "appl::grid::standard() using standard convolution" << std::endl; }

  CALCULATION  calculation() const { return m_type; }

  static std::string _calculation(CALCULATION C) { 
    switch (C) {
    case STANDARD:
      return "standard";
    case SHERPA:
      return "sherpa";
    case AMCATNLO:
      return "amcatnlo";
    case LAST_TYPE:
      return "last_type"; // NB: shouldn't ever be used
    }
    return "unknown";
  }

  /// reduce number of subprocesses if possible
  void shrink(const std::string& name, int ckmcharge=0);

  /// set bins to be combined after the convolution
  void combine( std::vector<int>& v) { if ( (m_combine=v).size() ) combineReference(true); }

  /// set combine the  be combined after the convolution
  void combineReference(bool force=false);

  void combineBins(std::vector<double>& v, int power=1 ) const;

  double fx(double x) const;
  double fy(double x) const;

  const appl_pdf* genpdf(int i) const { return m_genpdf[i]; }

  std::vector<double>&       userdata()       { return m_userdata; }
  const std::vector<double>& userdata() const { return m_userdata; }

protected:

  // internal common construct for the different types of constructor
  void construct(int Nobs,
		 int NQ2=50,  double Q2min=10000.0, double Q2max=25000000.0, int Q2order=4,  
		 int Nx=50,   double xmin=1e-5,     double xmax=0.9,         int xorder=3,
		 int order=2, 
		 std::string transform="f" );
  
protected:

  /// std::string manipulators to parse the pdf names 

  /// return chomped std::string
  static std::string chomptoken(std::string& s1, const std::string& s2)
  {
    std::string s3 = "";
    std::string::size_type pos = s1.find(s2);
    if ( pos != std::string::npos ) {
      s3 = s1.substr(0, pos);
      s1.erase(0, pos+1);
    }
    else { 
      s3 = s1.substr(0, s1.size());
      s1.erase(0, s1.size()+1);
    }
    return s3;
  } 
 
  static std::vector<std::string> parse(std::string s, const std::string& key) {
    std::vector<std::string> clauses;
    while ( s.size() ) clauses.push_back( chomptoken(s, key) );
    return clauses;
  }
  
  /// get the required pdf combinations from those registered   
  void findgenpdf( std::string s );

  /// add a generic pdf to the data base of registered pdfs
  void addpdf( const std::string& s, const std::vector<int>& combinations=std::vector<int>() );

  appl_pdf* genpdf(int i) { return m_genpdf[i]; }
  
public: 

  int subproc() const { return m_subproc; }

protected:

  // histograms for saving the observable
  TH1D*  m_obs_bins;
  TH1D*  m_obs_bins_combined;

  // order in alpha_s of tree level contribution 
  int  m_leading_order; 

  // how many orders in the calculation, lo, nlo, nnlo etc 
  int  m_order;

  // the actual weight grids themselves
  igrid** m_grids[MAXGRIDS]; /// up to MAXGRIDS grids LO, NLO, NNLO, Real virtual, etc 

  // total cross section qand uncertainty
  double m_total;
  double m_totalerror;

  // state variables
  double   m_run;
  bool     m_optimised;
  bool     m_trimmed;

  bool   m_normalised;

  bool   m_symmetrise; 
 
  // transform and pdf combination tags
  std::string m_transform; 
  std::string m_genpdfname; 

  // pdf combination class
  appl_pdf* m_genpdf[MAXGRIDS];

  static const std::string m_version;

  double m_cmsScale;

  double m_dynamicScale;

  /// bin by bin correction factors 
  //  std::vector<std::vector<double> > m_corrections;
  std::vector<correction>    m_corrections;
  std::vector<std::string>   m_correctionLabels;
  

  /// should we apply the corrections?
  bool m_applyCorrections;

  /// flag vector to determine whether each individual 
  /// correction should be applied
  std::vector<bool>  m_applyCorrection;

  std::string m_documentation;
  
  std::vector<double>                m_ckmsum;
  std::vector<std::vector<double> >  m_ckm2;
  std::vector<std::vector<double> >  m_ckm;

  CALCULATION     m_type; 

  bool            m_read;

  std::vector<int> m_combine;

  int  m_subproc;
  int  m_bin;

  std::vector<double> m_userdata;

};


};

// shouldn't have this, grid is too large a structure 
// grid operator*(const double& d, const appl::grid& g) { return g*d; }

std::ostream& operator<<(std::ostream& s, const appl::grid& mygrid);



#endif // __APPL_GRID_H 
