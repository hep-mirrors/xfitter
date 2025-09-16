/* emacs: this is -*- c++ -*- */
/**
 **   @file    smooth.h        
 **                   
 **   @author  sutt
 **   @date    Sat  5 Apr 2025 20:58:33 BST
 **
 **   $Id: smooth.h, v0.0   Sat  5 Apr 2025 20:58:33 BST sutt $
 **
 **   Copyright (C) 2025 sutt (sutt@cern.ch)    
 **
 **/

#ifndef  SMOOTH_H
#define  SMOOTH_H

#include <iostream>
#include <exception>

#include "appl_grid/appl_grid.h"

#include "TH1D.h"
#include "TPad.h"


typedef std::vector<std::vector<double> > covariance_t; 


inline std::ostream& operator<<( std::ostream& s, const std::vector<double>& v ) {
  //  for ( size_t i=0 ; i<v.size() ; i++ ) s << "\t" << v[i];
  for ( size_t i=0 ; i<v.size() ; i++ ) {
    double shi = v[i];
    if ( shi>9999 )    shi = 1000*long(shi*0.001); 
    else shi = 0.001*long(shi*1000);
    s << "\t" << shi;
  }
    
  return s;
}


inline std::ostream& operator<<( std::ostream& s, const std::vector<std::vector<double> > & v ) {
  for ( size_t i=0 ; i<v.size() ; i++ ) s << "\t" << v[i] << "\n";
  return s;
}



inline covariance_t operator*( const covariance_t& m0,  const covariance_t& m1 ) {
  /// no size checks here numbers of rows of first should match columns in 
  /// the second and vice versa - here we just treat everything as square,
  /// and don't even bother to check
  covariance_t m(m0.size(), std::vector<double>(m0.size(),0) );
  for ( size_t i=0 ; i<m0.size() ; i++ ) {
    for ( size_t j=0 ; j<m0.size() ; j++ ) {
      for ( size_t k=0 ; k<m0.size() ; k++ ) m[i][j] += m0[i][k]*m1[k][j];
    }
  }
  return m;
}


inline std::vector<double> operator*( const std::vector<double>& v0,  const std::vector<double>& v1 ) {
  if ( v0.size()!=v1.size() ) throw std::exception(); // ("vector size mismatch");
  std::vector<double>  v(v0.size(),0);
  for ( size_t i=0 ; i<v0.size() ; i++ ) v[i] += v0[i]*v1[i];
  return v;
}




std::vector<double>  operator*( const covariance_t& m0,  const std::vector<double>& v0 );


template<typename T>
inline std::vector<double>& operator+=( std::vector<T> v, const std::vector<T>& v1) {
  if ( v.size()!=v1.size() ) throw std::exception(); // ("vector size mismatch");
  for ( size_t i=0 ; i<v.size() ; i++ ) v[i] += v1[i];
  return v;
}


template<typename T>
inline std::vector<double> operator+( const std::vector<T>& v0, const std::vector<T>& v1) {
  if ( v0.size()!=v1.size() ) throw std::exception(); // ("vector size mismatch");
  std::vector<T> v(v0.size());
  for ( size_t i=0 ; i<v.size() ; i++ ) v[i] = v0[i] + v1[i];
  return v;
}



template<typename T>
inline std::vector<double> operator-( const std::vector<T>& v0, const std::vector<T>& v1) {
  if ( v0.size()!=v1.size() ) throw std::exception(); // ("vector size mismatch");
  std::vector<T> v(v0.size());
  for ( size_t i=0 ; i<v.size() ; i++ ) v[i] = v0[i] - v1[i];
  return v;
}


inline covariance_t operator+( const covariance_t& v0,  const covariance_t& v1 ) {
  if ( v0.size()!=v1.size() ) throw std::exception(); // ("vector size mismatch");
  covariance_t  v(v0.size(),std::vector<double>(v0.size(),0));
  for ( size_t i=0 ; i<v0.size() ; i++ ) v[i] = v0[i]+v1[i];
  return v;
}


inline covariance_t operator-( const covariance_t& v0,  const covariance_t& v1 ) {
  if ( v0.size()!=v1.size() ) throw std::exception(); // ("vector size mismatch");
  covariance_t  v(v0.size(),std::vector<double>(v0.size(),0));
  for ( size_t i=0 ; i<v0.size() ; i++ ) v[i] = v0[i]-v1[i];
  return v;
}


inline covariance_t identity( const std::vector<double>& v ) {
  covariance_t m(v.size(), std::vector<double>(v.size(),0) );
  for ( size_t i=0 ; i<v.size() ; i++ ) m[i][i] = v[i];
  return m;
}


inline covariance_t identity_inv( const std::vector<double>& v ) {
  covariance_t m(v.size(), std::vector<double>(v.size(),0) );
  for ( size_t i=0 ; i<v.size() ; i++ ) m[i][i] = 1/v[i];
  return m;
}


inline appl::TH1D fabs( const appl::TH1D& h ) {
  appl::TH1D h0(h);
  for ( size_t i=h0.size() ; i-- ; ) h0.y(i) = std::fabs(h0.y(i)); 
  return h0;
}


inline appl::TH1D fabs( const appl::TH1D* h ) { return fabs(*h); }


inline std::string basename( std::string s, const std::string& e="" ) { 
  if ( s.find("/")==std::string::npos ) return s;
  s = s.substr( s.find_last_of("/")+1, s.size() );
  if ( e=="" ) return s;
  if ( s.find(e)!=s.size()-e.size() ) return s;
  return s.substr( 0, s.size()-e.size() );
}



inline std::string dirname( std::string s ) { 
  if ( s.find("/")==std::string::npos ) return ".";
  s = s.substr( 0, s.find_last_of("/"));
  if ( s=="" ) return "/";
  return s;
}



std::vector<double> gaussian( int fullwidth, int n=1001, double scale=1 );

extern double fitscale;



class smooth {
   
public:

  smooth() { } 
  
  smooth( const appl::TH1D& h, const appl::TH1D& hr, const std::string& label="" ) : m_h(h), m_hr(hr) {

    //    std::cout << "smooth: label: " << label << std::endl;
    
    int n = h.size();

    int kn = 15;
    
    while ( kn>(n/2+1) ) {
      kn = n/2;
      if ( kn%2 == 0 ) kn += 1;
    }

    int start=0;

#if 0
    for ( size_t i=1 ; i<hr.size() ; i++ ) {
      if ( std::fabs(hr.y(i)/hr.y(i-1))>4 ) start = i;
      else break;
    }
#endif
    
    std::cout << "smooth: start: " << start << std::endl; 
    
    std::vector<double> w = gaussian( kn, 1001, fitscale );

    //    std::cout << w << std::endl;

    appl::TH1D hsmooth = h;

    hsmooth.clear();
    
    covariance_t  t( n, std::vector<double>(n,0) );
    covariance_t  t2( n, std::vector<double>(n,0) );
    covariance_t cov( n, std::vector<double>(n,0) );
    covariance_t cor( n, std::vector<double>(n,0) );
	
    {
      const double* p = &w[0]+kn/2;
      
      int lim = kn/2;

      if ( start>0 ) {
	for ( int i=0 ; i<start ; i++ ) t[i][i] = 1;
      }

      for ( size_t i=start ; i<t.size() ; i++ ) {	
	for ( int k=-lim ; k<=lim ; k++ ) {
	  
	  int j = i+k;
	  
#if 1
	  if ( j<start ) j=start;
	  if ( j>=int(t.size()) ) j=t.size()-1;
#else
	  if ( j<0 ) j=i;
	  if ( j>=int(t.size()) ) j=i;
#endif

	  t[i][j] += p[k];
	  
	}
      }
         
      for ( size_t i=0 ; i<t.size() ; i++ ) {
	for ( size_t j=0 ; j<t[i].size() ; j++ ) t[i][j]  *= hr.y(i)/hr.y(j);
      }
      
      for ( size_t i=0 ; i<t.size() ; i++ ) {
	for ( size_t j=0 ; j<t[i].size() ; j++ ) t2[i][j] = t[i][j]*t[i][j];
      }
           
    }

    //    std::cout << "\n\ntransform: " << t << "\n" << std::endl;
    
    m_hsmoothed = hsmooth;

    m_t = t;

    std::vector<double> ye2 = h.ye()*h.ye();

    std::vector<double> sm   = t*h.y();
    std::vector<double> sme2 = t2*ye2;

    //    std::vector<double> sme2 = ye2;
    
    for ( size_t i=ye2.size() ; i-- ; ) sme2[i] = std::sqrt(sme2[i]);

    m_hsmoothed.y()  = sm;

    m_hsmoothed.ye() = sme2;
    
    //    std::cout << "t: " << t << std::endl;     
    //    std::cout << "h: " << h << std::endl; 
    //    std::cout << "smoothed: " << m_hsmoothed << std::endl; 

    for ( size_t i=0 ; i<h.size() ; i++ ) { 
      for ( size_t j=0 ; j<h.size() ; j++ ) { 
	for ( size_t k=0 ; k<h.size() ; k++ ) cov[i][j] += t[i][k]*t[j][k]*h.ye(k)*h.ye(k);
      }
    }

    for ( size_t i=0 ; i<h.size() ; i++ ) { 
      for ( size_t j=0 ; j<h.size() ; j++ ) cor[i][j] = cov[i][j]/(h.ye(i)*h.ye(j));
    }

    m_k.clear();
    m_k.resize(m_h.size());
    
    for ( size_t i=0 ; i<m_h.size() ; i++ ) {
      m_k[i] = m_hsmoothed.y(i)/m_h.y(i);
    }
    
    m_cov = cov;

#if 0
    for ( size_t i=0 ; i<cov.size() ; i++ ) {
      std::cout << "\t" << i
		<< "\t: " << m_hsmoothed.y(i)
		<< "\t"   << m_hsmoothed.ye(i)
		<< "\t"   << std::sqrt(cov[i][i])
		<< "\t (" << (m_hsmoothed.ye(i)-std::sqrt(cov[i][i]))/m_hsmoothed.ye(i) << ")" << std::endl; 
      
    }
#endif

    
  }


  smooth( const smooth& s) :
    m_h(s.m_h), m_hr(s.m_hr), m_hsmoothed(s.m_hsmoothed),
    m_k(s.m_k), m_cov(s.m_cov) { }  

  
  virtual ~smooth() { } 

  operator appl::TH1D() const { return m_hsmoothed; }
  
  const std::vector<double> k() const { return m_k; }

  covariance_t cov() const { return m_cov; }

  covariance_t transform() const { return m_t; }

private:

  appl::TH1D m_h;
  appl::TH1D m_hr;

  appl::TH1D m_hsmoothed;

  std::vector<double> m_k;

  covariance_t m_cov;

  covariance_t m_t;

};

inline std::ostream& operator<<( std::ostream& s, const smooth& _ ) { 
  return s;
}



#endif  // SMOOTH_H 










