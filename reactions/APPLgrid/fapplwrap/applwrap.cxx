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


#include "fit.h"


bool applwrap::m_fastsmooth = true;
int  applwrap::m_ratiobase = 0;



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
      if ( std::fabs(c[i][j])<1e-20 ) std::printf("      ---");
      else                            printf("   %6lf", c[i][j] );
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
  
  //  std::vector<std::vector<double>  > transform = mg->getCovariance();
  covariance_t transform = mg->getCovariance();

  if ( nloops!=mg->nloops() )  return mg->aconvolute( pdf1, pdf2, alphas, nloops, rscale, fscale, escale );

  if ( transform.size() == 0 && m_fitorder==-1 ) return mg->aconvolute( pdf1, pdf2, alphas, nloops, rscale, fscale, escale );

  /// check that we have order-by-order references ...

  for ( int iloops=mg->nloops()+1 ; iloops-- ; ) {
      xr[iloops] = mg->getReference(iloops);
      if ( xr[iloops]->size()==0 ) return mg->aconvolute( pdf1, pdf2, alphas, nloops, rscale, fscale, escale );
  }
  
  /// perform the convolutions order by order

  for ( int iloops=mg->nloops()+1 ; iloops-- ; ) {
    
    xs[iloops] = mg->aconvolute( pdf1, pdf2, alphas, -1*iloops, rscale, fscale, escale );

    std::vector<double> ye = xs[iloops]->y()*xr[iloops]->ye()/xr[iloops]->y();

    xs[iloops]->ye() = ye;
    
  }
  
  /// smooth the nnlo component with the kernel stored in the
  /// covariance matrix


  //  appl::TH1D xsum = *xs[0];

  //  for ( int i=1 ; i<mg->nloops()+1 ; i++ ) xsum += *xs[i];

  appl::TH1D* sxsum = new appl::TH1D();
    
  if ( m_fitorder>-1 ) { 

    /// fit with a polynomial ...

    /// FIXME: need to put in some limit to the fitted region to avoid
    ///        discontinuities as in first bin of the dijet grids  
    
#if 0
    std::cout << "-------------------------------------" << std::endl;
    std::cout << " perform polynomial fit " << std::endl;
    std::cout << "-------------------------------------" << std::endl;
#endif
  
    //  appl::TH1D hratio1 = *xs[1]/xr[0]->y();
    appl::TH1D hratio2 = *xs[2]/xr[0]->y();
    
    //  appl::TH1D cs1 = generate( hratio1, m_fit_order ) ; 
    appl::TH1D cs2 = generate( hratio2, m_fitorder ) ; 
    
    //  cs1 *= xr[0]->y();
    cs2 *= xr[0]->y();
    
    *sxsum = cs2;  /// smoothed NNLO component
    
    *sxsum += *xs[0] + *xs[1]; /// add the LO and NLO
    
    /// calculate the variance ...

    for ( size_t i=0 ; i<xs[0]->size() ; i++ ) {
      sxsum->ye(i) =  std::sqrt( xs[0]->ye(i)*xs[0]->ye(i) + xs[1]->ye(i)*xs[1]->ye(i) + cs2.ye(i)*cs2.ye(i) );
    }
  }
  else {  

    /// smoothing ....

#if 0
    std::cout << "-------------------------------------" << std::endl;
    std::cout << " perform GK smoothing " << std::endl;
    std::cout << "-------------------------------------" << std::endl;
#endif
    
    appl::TH1D cs = *xs[2];

    /// covariance ???
    
    covariance_t   scov( cs.size(), std::vector<double>( cs.size(), 0 ) );
    
    if ( m_fastsmooth ) { 
      
      /// use the smoothing trasform encoded in the grid ...
      
      /// smooth the NNLO component ...
      
      cs.y() = transform*xs[2]->y();
      
      /// calculate the covariance ...
      
      for ( size_t i=0 ; i<cs.size() ; i++ ) { 
	for ( size_t j=0 ; j<cs.size() ; j++ ) { 
	  for ( size_t k=0 ; k<cs.size() ; k++ ) scov[i][j] += transform[i][k]*transform[j][k]*cs.ye(k)*cs.ye(k);
	}
      }
      
    }
    else { 
      
      
      smooth sm;
      
      /// use just the LO component ...
      /// generate the smoothing trasform and smooth the grid ...
      if      ( ratiobase()==0 ) sm = smooth( *xs[2],  *xs[0] );
      else if ( ratiobase()==1 ) sm = smooth( *xs[2], (*xs[0])+(*xs[1]) ); 
      else std::cerr << "Incorrect smoothing base ratio" << std::endl;
      
      /// get the smoothed output ...
      
      cs = sm;
      
      /// and the covariance
      
      scov = sm.cov();
      
      /// store in the grid if required ...
      
      covariance_t trans = sm.transform();
      
      mg->getCovariance() = trans;
      
    }
    
    *sxsum = cs;  /// smoothed NNLO component
    
    *sxsum += *xs[0] + *xs[1]; /// add the LO and NLO
       
    for ( size_t i=0 ; i<xs[0]->size() ; i++ ) {
      sxsum->ye(i) =  std::sqrt( xs[0]->ye(i)*xs[0]->ye(i) + xs[1]->ye(i)*xs[1]->ye(i) + cs.ye(i)*cs.ye(i) );
    }

    /// calculate the full covariance ...
    
    for ( size_t i=scov.size() ; i-- ; ) {
      scov[i][i] += xs[0]->ye(i)*xs[0]->ye(i) + xs[1]->ye(i)*xs[1]->ye(i);
      sxsum->ye(i) = std::sqrt(scov[i][i]);
    }
    
    m_scovariance = scov;
  }
   
  //  print_covariance(scov);
  
  //  print_covariance(scov, sxsum->ye());
  
  /// add the LO and NLO diagonal components ...

  return sxsum; 
}


