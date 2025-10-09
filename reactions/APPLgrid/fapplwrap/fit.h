/* emacs: this is -*- c++ -*- */
/**
 **   @file    fit.h        
 **                   
 **   @author  sutt
 **   @date    Fri 3 Oct 2025 12:12:44 BST
 **
 **   $Id: fit.h, v0.0   Fri 3 Oct 2025 12:12:44 BST sutt $
 **
 **   Copyright (C) 2025 sutt (sutt@cern.ch)    
 **
 **/



#ifndef  APPLWRAP_FIT_H
#define  APPLWRAP_FIT_H

#include "appl_grid/appl_TH1D.h"

// #include "simpletimer.h"

#include "solve.h"


#ifdef USE_EIGEN 
#include "Eigen/Dense"
using namespace Eigen;
#endif

double var( double x, const matrix_t& cov, const vector_t& a ) {

  vector_t phi(a.size(),0);

  double x0 = 1;

  //  std::cout << "a.size(): " << a.size() << std::endl;
  
  for ( size_t i=0 ; i<=a.size() ; i++ ) { 

    phi[i] = x0;
    
    x0 *= x;
    
  }
    
  double var = dot( phi, (cov*phi) );

  return var;

}


vector_t pol( double x, const matrix_t& cov, const vector_t a ) {

  vector_t phi(a.size(),0);
  
  double x0 = 1;
  
  for ( size_t i=0 ; i<a.size() ; i++ ) { 

    phi[i] = x0;
    
    x0 *= x;
    
  }

  vector_t f(2,0);

  f[0] = dot(phi,a);
  
  f[1] = sqrt( dot( phi, (cov*phi) ) );

  return f;

}


void symmetrise( matrix_t& m ) {
  for ( size_t i=m.size() ; i-- ; ) { 
    for ( size_t j=i ; j-- ; ) {
      m[i][j] = 0.5*(m[i][j] + m[j][i]);
      m[j][i] = m[i][j];
    }
  }
}



double sigma( double x, const matrix_t& cov, const vector_t a ) {
  return std::sqrt(var( x, cov, a ));
}


#ifdef USE_EIGEN 

MatrixXd convert( const matrix_t& m ) {
  MatrixXd me(m.size(), m.size());
  for ( size_t i=0 ; i<m.size() ; i++ ) { 
    for ( size_t j=0 ; j<m.size() ; j++ ) me(i,j) = m[i][j];
  }
  return me;
}

matrix_t convert( const MatrixXd& me ) {
  matrix_t m(me.rows(), vector_t(me.cols()));
  for ( size_t i=me.rows() ; i-- ;  ) { 
    for ( size_t j=me.cols() ; j-- ; ) m[i][j] = me(i,j); 
  }
  return m;
}

vector_t convert( const VectorXd& ve ) {
  vector_t v(ve.rows());
  for ( size_t i=ve.rows() ; i-- ; ) v[i] = ve(i); 
  return v;
}

VectorXd  convert( const vector_t& v ) {
  VectorXd  ve(v.size());
  for ( size_t i=v.size() ; i-- ; ) ve(i) = v[i]; 
  return ve;
}

#endif


appl::TH1D generate( const appl::TH1D& hin, size_t N ) {

  //  std::cout << "\n-------------------\ngenerate() order N: " << N << std::endl;
  
  vector_t x(hin.size());
  for ( size_t i=0 ; i<hin.size() ; i++ ) x[i] = std::log10(0.5*(hin.lo(i)+hin.hi(i)));

  //  std::cout << "x:  " << x << std::endl;
  
  vector_t y  = hin.y();
  vector_t ye = hin.ye();
  
  matrix_t A(N+1, vector_t(N+1,0) );

  vector_t s( 2*N+1, 0 );
  vector_t c(   N+1, 0 );

  vector_t w(ye.size());

  for ( size_t i=0 ; i<w.size() ; i++ ) w[i] = 1/(ye[i]*ye[i]);
  
  for ( size_t k=0 ; k<x.size() ; k++ ) {

    double s0 = w[k];
      
    for ( size_t i=0 ; i<=2*N ; i++ ) { 
      s[i] += s0;
      s0   *= x[k];
    }
    
    double y0 = w[k]*y[k];

    for ( size_t i=0 ; i<=N ; i++ ) {
      c[i] += y0;
      y0 *= x[k];
    }
    
  }
      
  
  for ( size_t i=0 ; i<=N ; i++ ) {
    A[i][i] = s[2*i];
    for ( size_t j=i+1 ; j<=N ; j++ ) {
      A[i][j] = s[i+j];
      A[j][i] = s[i+j];
    }
  }

#ifndef USE_EIGEN
  
  solve ss( A, c );

  vector_t soln = ss.solution();
  matrix_t inv  = ss.inverse(); 

#else 

  //  std::cerr << "USING EIGEN" << std::endl;
  
  MatrixXd AM = convert( A );
  VectorXd cM = convert(c);

  ///  MatrixXd AMinv = AM.inverse();
  ///  matrix_t inv = convert( AMinv );
  
  /// solving using the inverse is not efficient, but this was
  /// just for a check wiht eigen
  /// vector_t soln = inv*c;
  
  Eigen::LDLT<MatrixXd> ldlt(AM);

  VectorXd sM = ldlt.solve(cM);
  
  MatrixXd AMinv = ldlt.solve(MatrixXd::Identity(AM.rows(), AM.cols()));

  vector_t soln = convert(sM);
  matrix_t inv  = convert(AMinv);
  
#endif

  //  std::cout << "inv:  " << inv << std::endl;
  //  std::cout << "soln: " << soln << std::endl;
  
  size_t Nsize = x.size();
  
  vector_t  yy(Nsize,0);
  vector_t yye(Nsize,0);


  appl::TH1D h = hin;
  
  for ( size_t i=0 ; i<x.size() ; i++ ) {

    double xx = x[i];

    vector_t f = pol( xx, inv, soln );

    yy[i]  = f[0];
    yye[i] = f[1];
  }
  
  h.y()  = yy;
  h.ye() = yye;

  return h;
}






#endif // APPLWRAP_FIT_H

