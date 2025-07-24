/**
 **     @file    fapplwrap.cxx
 **
 **     @brief   fortran callable wrapper functions for the c++  
 **              appl grid project applwrap smoother
 **
 **     @author  mark sutton
 **     @date    Wed May 21 14:31:36 CEST 2008 
 **
 **     @copyright (C) 2002 mark sutton (sutt @ cern.ch) 
 **
 **     $Id: fappl.cxx, v1.0   Wed May 21 14:31:36 CEST 2008 sutt $
 **
 **/

/// externally defined alpha_s and pdf routines for fortran 
/// callable convolution wrapper
extern "C" double fnalphas_(const double& Q); 
extern "C" void   fnpdf_(const double& x, const double& Q, double* f);

#include "applwrap.h" 


/// externally defined alpha_s and pdf routines for fortran 
/// callable convolution wrapper
extern "C" double fnalphas_(const double& Q); 
extern "C" void   fnpdf_(const double& x, const double& Q, double* f);



extern std::map<int,appl::grid*> _grid;

void throw_exception( const std::string& msg, int id, const std::string& s="" );




extern "C" void sconvolutewrap_(const int& id, double* data, double* dataerr, 
				void (*pdf)(const double& , const double&, double* ),  
				double (*alphas)(const double& ) ) {  

  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) { 
    appl::grid*    g = gitr->second;
    applwrap ag(g);
    appl::TH1D* v = ag.sconvolute( pdf, alphas );
    for ( unsigned i=0 ; i<v->size() ; i++ ) {
      data[i]    = v->y(i);
      dataerr[i] = v->ye(i);
    }
    delete v;
  }
  else throw_exception( "No grid with id ", id );  
}


extern "C" void sconvolute_(const int& id, double* data, double* dataerr) {
  sconvolutewrap_(id, data, dataerr, fnpdf_, fnalphas_); 
}



extern "C" void sfullconvolutewrap_(const int& id, double* data, double* dataerr, 
				    void (*pdf)(const double& , const double&, double* ),  
				    double (*alphas)(const double& ),
				    const int& nloops,
				    const double& rscale, const double& fscale  ) {  

  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) { 
    appl::grid*    g = gitr->second;
    appl::TH1D*    v = applwrap(g).sconvolute( pdf, alphas, nloops, rscale, fscale); 
    for ( unsigned i=0 ; i<v->size() ; i++ ) {
      data[i]    = v->y(i);
      dataerr[i] = v->ye(i);
    }
    delete v;
  }
  else throw_exception( "No grid with id ", id );  
}


extern "C" void sfullconvolute_(const int& id, double* data, double* dataerr, 
				const int& nloops,
				const double& rscale, const double& fscale  ) {  
  sfullconvolutewrap_( id, data, dataerr, fnpdf_, fnalphas_, nloops, rscale, fscale);
}




extern "C" void sconvolute_covariancewrap_(const int& id, double* data, double* dataerr, 
				void (*pdf)(const double& , const double&, double* ),  
				double (*alphas)(const double& ) ) {  

  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) { 
    appl::grid*    g = gitr->second;
    applwrap ag(g);
    appl::TH1D* v = ag.sconvolute( pdf, alphas );
    for ( unsigned i=0 ; i<v->size() ; i++ ) {
      data[i]    = v->y(i);
      dataerr[i] = v->ye(i);
    }
    delete v;
  }
  else throw_exception( "No grid with id ", id );  
}


extern "C" void sconvolute_covariance_(const int& id, double* data, double* dataerr) {
  sconvolute_covariancewrap_(id, data, dataerr, fnpdf_, fnalphas_); 
}

