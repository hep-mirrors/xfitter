/* emacs: this is -*- c++ -*- */
/**
 **   @file    applwrap.h        
 **                   
 **   @author  sutt
 **   @date    Sun 13 Jul 2025 20:17:37 BST
 **
 **   $Id: applwrap.h, v0.0   Sun 13 Jul 2025 20:17:37 BST sutt $
 **
 **   Copyright (C) 2025 sutt (sutt@cern.ch)    
 **
 **/


#ifndef  APPLWRAP_H
#define  APPLWRAP_H

#include <iostream>

#include "appl_grid/appl_grid.h"
#include "appl_grid/histogram.h"


typedef std::vector<std::vector<double> > covariance_t; 




class applwrap { // : public appl::grid {

public:


  /// better make sure that the passed in grid is not deleted ...
  applwrap( appl::grid* g ) : mmanage(false), mg(g), m_fitorder(-1) { }

  /// better make sure that the passed in grid is not deleted ...
  applwrap( appl::grid& g ) : mmanage(false), mg(&g), m_fitorder(-1) { }

  applwrap( const std::string& s ) : mmanage(true), m_fitorder(-1) {
    std::cout << "applwrap: reading new grid: " << s << std::endl; 
    mg = new appl::grid(s);
  }
    
  virtual ~applwrap() { if ( mmanage ) delete mg; } 


  /// control fitting or smoothing ...
  /// default   order = -1, run smoothing
  /// otherwise order > -1 fit polynomial, 0 is a constant, 1 is linear
  int  fit_order() const  { return m_fitorder; }
  void fit_order( int i ) { m_fitorder=i; }


  /// run the "smoothed" or "fitted" convolution ...
  
  appl::TH1D* sconvolute(void   (*pdf)(const double& , const double&, double* ), double (*alphas)(const double& ) ) { 
    return sconvolute( pdf, alphas, mg->nloops() ); 
  }
  
  appl::TH1D* sconvolute(void (*pdf)(const double& , const double&, double* ), double (*alphas)(const double& ),
			 const int& nloops,
			 const double& rscale=1, const double& fscale=1, double escale=1 ) {
    return sconvolute( pdf, pdf, alphas, nloops, rscale, fscale, escale );
  }
  
  appl::TH1D* sconvolute(void   (*pdf1)(const double& , const double&, double* ),
			 void   (*pdf2)(const double& , const double&, double* ),
			 double (*alphas)(const double& ),
			 const int& nloops,
			 const double& rscale=1, const double& fscale=1, double escale=1 ); 
  

  appl::TH1D* aconvolute(void (*pdf)(const double& , const double&, double* ), double (*alphas)(const double& ) ) { 
    appl::TH1D* h  = mg->aconvolute( pdf, alphas );
    appl::TH1D* hr = mg->getReference();
    for ( size_t i=0 ; i<h->size() ; i++ ) h->ye(i) = hr->ye(i)*h->y(i)/hr->y(i);
    return h;
  }

  std::vector<std::vector<double> >& getCovariance() { return mg->getCovariance(); }
  
  appl::grid* g() { return mg; }

  covariance_t covariance() const { return m_scovariance; }

 
protected:  

  static bool& fastsmooth() { return m_fastsmooth; }
  static void  fastsmooth(bool b) { m_fastsmooth = b; }

  static int   ratiobase()        { return m_ratiobase; }
  static void  ratiobase( int i ) { m_ratiobase=i; }
 
private:

  bool mmanage;

  appl::grid* mg;

  covariance_t m_scovariance;

  int          m_fitorder;  
  
  static bool  m_fastsmooth;
  static int   m_ratiobase;

};


inline std::ostream& operator<<( std::ostream& s, const applwrap& _a ) { 
  return s;
}


void print_covariance( const covariance_t&  c, const std::vector<double>& v);

void print_covariance( const covariance_t&  c );


#endif  // APPLWRAP_H 










