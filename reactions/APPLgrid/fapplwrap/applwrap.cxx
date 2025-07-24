/**
 **   @file    applwrap.cxx         
 **   
 **
 **   @author sutt
 **   @date   Sun 13 Jul 2025 20:19:39 BST
 **
 **   $Id: applwrap.cxx, v0.0   Sun 13 Jul 2025 20:19:39 BST sutt $
 **
 **   Copyright (C) 2025 sutt (sutt@cern.ch)    
 **/

#include "applwrap.h"
#include "appl_grid/histogram.h"

#include "smooth.h"


template<typename T>
std::vector<T> operator*( const std::vector<T>& v0, const std::vector<T>& v1 ) {

  if ( v0.size()!=v1.size() ) throw std::runtime_error("vectors do not match");
  
  std::vector<T> v(v0);

  for ( size_t i=v0.size() ; i-- ; ) v[i] *=v1[i];

  return v;
}


template<typename T>
std::vector<T> operator/( const std::vector<T>& v0, const std::vector<T>& v1 ) {

  if ( v0.size()!=v1.size() ) throw std::runtime_error("vectors do not match");
  
  std::vector<T> v(v0);

  for ( size_t i=v0.size() ; i-- ; ) if ( v1[i]!=0 ) v[i] /= v1[i]; else throw std::runtime_error("division by 0");

  return v;
}


void print_covariance( const covariance_t&  c, const std::vector<double>& v) {
  for ( size_t i=0 ; i<10 && i<c.size() ; i++ ) { 
    for ( size_t j=0 ; j<10 && j<c.size() ; j++ ) { 
      if ( std::fabs(c[i][j])<1e-20 ) printf("      ---");
      else                            printf("  %7.5lf",c[i][j]/(v[i]*v[j]));
    }
    if ( c.size()>10 ) printf("     ...");
    printf("\n");
  }
  if ( c.size()>10 ) printf("     ...");
  printf("\n");
}


void print_covariance( const covariance_t&  c ) {
  for ( size_t i=0 ; i<10 && i<c.size() ; i++ ) { 
    for ( size_t j=0 ; j<10 && j<c.size() ; j++ ) { 
      if ( std::fabs(c[i][j])<1e-20 ) printf("     ---");
      else                            printf("  %6lf",c[i][j]);
    }
    if ( c.size()>10 ) printf("     ...");
    printf("\n");
  }
  if ( c.size()>10 ) printf("     ...");
  printf("\n");
}



  
appl::TH1D* applwrap::sconvolute(void (*pdf1)(const double& , const double&, double* ), 
				 void (*pdf2)(const double& , const double&, double* ),
				 double (*alphas)(const double& ),
				 const int& nloops,
				 const double& rscale, const double& fscale, double escale ) {
  
  if ( mg==0 ) std::cerr << "applwrap whoops: grid not setup" << std::endl; 
  
  appl::TH1D* xs[3];  /// convolutions order by order
  appl::TH1D* xr[3];  /// references order by order
  
  std::vector<std::vector<double>  > cov = mg->getCovariance();

  if ( nloops!=mg->nloops() ) return mg->aconvolute( pdf1, pdf2, alphas, nloops, rscale, fscale, escale );
  
  if ( cov.size() == 0 )    return mg->aconvolute( pdf1, pdf2, alphas, nloops, rscale, fscale, escale );

  /// perform the convolutions order by order
  
  for ( int iloops=mg->nloops()+1 ; iloops-- ; ) {

    xs[iloops] = mg->aconvolute( pdf1, pdf2, alphas, -1*iloops, rscale, fscale, escale );

    xr[iloops] = mg->getReference(iloops);
    
    std::vector<double> ye = xs[iloops]->y()*xr[iloops]->ye()/xr[iloops]->y();
      
    xs[iloops]->ye() = ye;
    
  }
  
  /// smooth the nnlo component with the kernel stored in the
  /// covariance matrix
  
  appl::TH1D xsum = *xs[0] + *xs[1] + *xs[2];

  smooth sm( *xs[2], *xs[0] );

  appl::TH1D* sxsum = new appl::TH1D();
  
  *sxsum = sm;  /// smoothed NNLO component
  
  *sxsum += *xs[0] + *xs[1]; /// add the LO and NLO

  /// cobvariance ???

#if 0  
  std::cout << *sxsum << std::endl;

  for ( size_t i=sxsum->size() ; i-- ; ) {
    std::cout << "sxsum2 " << i << "\t" << (sxsum->y(i)*sxsum->y(i)) << std::endl;
  }
#endif
  
  m_covariance = sm.cov();

  covariance_t& scov = m_covariance;

  //  print_covariance(scov);
  
  //  print_covariance(scov, sxsum->ye());

  for ( size_t i=scov.size() ; i-- ; ) {
    scov[i][i] += xs[0]->ye(i)*xs[0]->ye(i) + xs[1]->ye(i)*xs[1]->ye(i);
  }

  //  print_covariance(scov, sxsum->ye());

  
#if 0
  appl::TH1D* ref = mg->getReference();

  /// these are the scaled uncertainties for the nominal convolution
  /// do not need to be computed her as they are never actually used
   
  for ( size_t i=0 ; i<ref->size() ; i++ ) {

    double er = ref->ye(i)/ref->y(i);

    xsum.ye(i) = er * xsum.y(i);
    
  }  
#endif
  
  return sxsum; 
}


