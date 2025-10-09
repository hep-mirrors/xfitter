/* emacs: this is -*- c++ -*- */
/**
 **   @file    solve.h        
 **                   
 **   @author  sutt
 **   @date    Fri 26 Sep 2025 12:12:44 BST
 **
 **   $Id: solve.h, v0.0   Fri 26 Sep 2025 12:12:44 BST sutt $
 **
 **   Copyright (C) 2025 sutt (sutt@cern.ch)    
 **
 **/


#ifndef  APPLWRAP_SOLVE_H
#define  APPLWRAP_SOLVE_H

#include "amconfig.h"

#include <iostream>
#include <vector>

typedef std::vector<std::vector<double> >  matrix_t; 
typedef std::vector<double>  vector_t; 

#if 0
std::ostream& operator<<( std::ostream& s, const vector_t& v ) { 
  for ( size_t i=0 ; i<v.size() ; i++ ) s << "\t" << v[i];
  return s;
}

std::ostream& operator<<( std::ostream& s, const matrix_t& m ) { 
  for ( size_t i=0 ; i<m.size() ; i++ ) s << m[i] << "\n";
  return s;
}
#endif

inline vector_t& operator*=( vector_t& v, const double& d ) {
  for ( size_t i=v.size() ; i-- ; ) v[i] *= d;
  return v;
}



inline vector_t operator*( vector_t v, const double& d ) {
  for ( size_t i=v.size() ; i-- ; ) v[i] *= d;
  return v;
}


inline vector_t operator*( const double& d, vector_t v ) {
  for ( size_t i=v.size() ; i-- ; ) v[i] *= d;
  return v;
}


inline vector_t& operator-=( vector_t& v, const vector_t& v0 ) {
  for ( size_t i=v.size() ; i-- ; ) v[i] -= v0[i];
  return v;
}

inline vector_t operator-( vector_t v, const vector_t& v0 ) {
  for ( size_t i=v.size() ; i-- ; ) v[i] -= v0[i];
  return v;
}


inline vector_t& operator+=( vector_t& v, const vector_t& v0 ) {
  for ( size_t i=v.size() ; i-- ; ) v[i] += v0[i];
  return v;
}

inline vector_t operator+( vector_t v, const vector_t& v0 ) {
  for ( size_t i=v.size() ; i-- ; ) v[i] += v0[i];
  return v;
}


inline vector_t& operator*=( vector_t& v, const vector_t& v0 ) {
  for ( size_t i=v.size() ; i-- ; ) v[i] *= v0[i];
  return v;
}

inline vector_t& operator/( vector_t& v, const vector_t& v0 ) {
  for ( size_t i=v.size() ; i-- ; ) v[i] /= v0[i];
  return v;
}


/// this is what it should really be ...
// double operator*( const vector_t& v, const vector_t& v0 ) {
//  double s = 0;
//  for ( size_t i=v.size() ; i-- ; ) s += v[i]*v0[i];
//  return s;
// }

double dot( const vector_t& v, const vector_t& v0 ) {
  double s = 0;
  for ( size_t i=v.size() ; i-- ; ) s += v[i]*v0[i];
  return s;
}



double sum( const vector_t& v ) {
  double s = 0;
  for ( size_t i=v.size() ; i-- ; ) s += v[i];
  return s;
}


void print( const matrix_t& m, const std::string& s="" ) {
  if ( s!="" ) std::cout << s << std::endl;
  for ( size_t i=0 ; i<m.size() ; i++ ) { 
    for ( size_t j=0 ; j<m[i].size() ; j++ ) {
      if ( j==m.size() ) std::cout << "\t";
      if ( j>15 ) { std::cout << "\t..."; break; }
      std::cout << "\t" << m[i][j];
    }
    std::cout << std::endl; 
  }
}


inline matrix_t operator*( matrix_t m, const double& d ) {
  for ( size_t i=m.size() ; i-- ; ) m[i] *= d;
  return m;
}

inline matrix_t operator*( const double d, matrix_t m ) {
  for ( size_t i=m.size() ; i-- ; ) m[i] *= d;
  return m;
}


inline matrix_t operator*=( matrix_t& m, const double& d ) {
  for ( size_t i=m.size() ; i-- ; ) m[i] *= d;
  return m;
}



#if 0
std::ostream& operator<<( std::ostream& s, const matrix_t& m ) {
  for ( size_t i=0 ; i<m.size() ; i++ ) { 
    for ( size_t j=0 ; j<m[i].size() ; j++ ) {
      if ( j>15 ) { s << "\t..."; break; }
      s << "\t" << m[i][j];
    }
    s << "\n"; 
  }
  return s;
}

std::ostream& operator<<( std::ostream& s, const vector_t& v ) {
  for ( size_t i=0 ; i<v.size() ; i++ ) {
    s << "\t" << v[i] << "\n";
  }
  return s;
}
#endif


inline vector_t operator*( const matrix_t& m, const vector_t& v ) {
  vector_t vo(v.size(),0);
  for ( size_t i=v.size() ; i-- ; ) vo[i] += dot(m[i],v);
  return vo;
}


#if 0
inline matrix_t operator*( const matrix_t& m, const matrix_t& m0 ) {
  matrix_t mo(m.size(),vector_t(m.size(),0));
  vector_t vo(m.size(),0);
  for ( size_t i=m.size() ; i-- ; ) { 
    for ( size_t j=m.size() ; j-- ; ) {
      for ( size_t k=m.size() ; k-- ; ) mo[i][j] += m[i][k]*m0[k][j];
    }
  }
  return mo;
}
#endif



class solve {

public:

  solve( matrix_t m, const vector_t& v ) : m_m(m), m_v(v) {

    size_t n = m.size();
    
    for ( size_t i=0; i<m.size() ; i++ ) {
      m[i].push_back( v[i] );
      for ( size_t j=0 ; j<n ; j++ ) m[i].push_back(0);
      m[i][n+1+i] = 1;
    }

    /// gauss-jordon elimination
    
    for ( size_t i=m.size() ; i-- ; ) {
      m[i] *= 1/m[i][i];

      for ( size_t j=i ; j-- ;  ) {

	double c = m[j][i]; // /m[i][i];

	m[j] = (m[j]-(c*m[i]));

	m[j][i] = 0; /// not really needed but just in case
      }
      
    }

    /// forward substitution
    
    //    for ( size_t i=m.size() ; i-- ; ) {
    for ( size_t i=0 ; i<m.size() ; i++ ) {

      for ( size_t j=i+1 ; j<m.size() ; j++ ) { 
      
	double c = m[j][i]; // /m[i][i];
      
	m[j] -= (c*m[i]);
	
	m[j][i] = 0; /// not really needed but just in case
      }

    }

    /// since I do this before the gauss-jordan elimination,
    /// I avoid the need to do any subsequent divisions
    // for ( size_t i=0 ; i<m.size() ; i++ ) m[i] *= 1/m[i][i];
    
    //    print(m, "\nback substuted");
    
    vector_t soln(m.size());

    for ( size_t i=m.size() ; i-- ; ) soln[i] = m[i][m.size()];

    m_y = soln;
    

    /// store the inverse as it will be needed later ...
    
    m_inv = matrix_t(m_m.size());

    for ( size_t i=0 ; i<m_m.size() ; i++ ) {
      m_inv[i].resize(m_m.size());
      for ( size_t j=0 ; j<m_m.size() ; j++ ) {
	m_inv[i][j] = m[i][m_m.size()+1+j];
      }
    }
    
  } 

  
  virtual ~solve() { } 

  vector_t solution() const { return m_y; }
  matrix_t inverse()  const { return m_inv; }
  
private:

  matrix_t m_m;
  vector_t m_v;

  vector_t m_y;

  matrix_t m_inv;

};


inline std::ostream& operator<<( std::ostream& s, const solve& _s ) { 
  return s;
}


#endif  // APPLWRAP_SOLVE_H 










