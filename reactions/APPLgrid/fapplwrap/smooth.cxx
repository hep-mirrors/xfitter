/**
 **   @file    smooth.cxx         
 **   
 **
 **   @author sutt
 **   @date   Sat  5 Apr 2025 20:58:24 BST
 **
 **   $Id: smooth.cxx, v0.0   Sat  5 Apr 2025 20:58:24 BST sutt $
 **
 **   Copyright (C) 2025 sutt (sutt@cern.ch)    
 **/


#include <cmath>
#include <sys/stat.h>

#include "appl_grid/appl_file.h"


#include "smooth.h"

appl::TH1D*  smooth::REF = 0;
covariance_t smooth::T   = covariance_t(); 
covariance_t smooth::T2  = covariance_t(); 





appl::TH1D* get_histogram( const std::string& filename ) {

  struct stat _fileinfo;

  appl::TH1D* ref = 0;
  
  if ( stat(filename.c_str(),&_fileinfo) )   {    
    throw appl::grid::exception( std::string("grid::grid() cannot find file ") + filename  ); 
  }

  // std::cout << "reading grid from file " << filename;
  
  appl::file gridfile( filename, "r" );

  if ( !gridfile.isopen() ) {
    // std::cout << std::endl;
    throw appl::grid::exception( std::string("grid::grid() cannot open file: ") + filename );
  }


  if ( gridfile.index().find( "reference_internal" ).size > 0 ) { 
    ref = new appl::TH1D( gridfile.Read<appl::TH1D>( "reference_internal" ) );
  }
  else { 
    if ( gridfile.index().find( "reference" ).size > 0 ) { 
      ref = new appl::TH1D( gridfile.Read<appl::TH1D>( "reference" ) );
    }
    else { 
      throw appl::grid::exception( "cannot read reference histogram" );
    }
  }
  
  return ref;
}






const double inorm = 1/std::sqrt(2*M_PI);

double gauss( double x, double sigma ) {
  double norm = inorm/sigma;
  return norm*std::exp( -x*x/(2*sigma*sigma) );
}



std::vector<double> gaussian( int fullwidth, int n, double scale ) {

  if ( fullwidth==1 ) return std::vector<double>(1,1);
  
  int N = 2*(fullwidth/2)+1;

  if ( 2*(n/2)==n ) n += 1;
  
  int Nh = N*n;

  std::vector<double> h(Nh,0);

  double imid = Nh/2; 
  
  double sigma = 0.1*fullwidth;

  double integral0 = 0;

  //  std::cout << "SUTT gaussian: N: " << Nh << std::endl;
  
  for ( size_t i=0 ; i<h.size() ; i++ ) {

    double x = double(i)-imid;
    
    h[i] = gauss( x, sigma*n*scale);

    integral0 += h[i];

    //    std::cout << i << "  " << x << "\t" << h[i] << std::endl;

  }

  double diff = 0.5*(1-integral0);
  
  //  std::cout << "SUTT integral: " << integral0 << std::endl;
  
  std::vector<double> w(N,0);

  w[0] += diff;
  w[w.size()-1] += diff;
  
  double integral = 0;
  
  for ( size_t i=0 ; i<w.size() ; i++ ) {

    for ( int j=0 ; j<n ; j++ ) w[i] += h[i*n+j];
    
    //   std::cout << i << "  " << "\t" << w[i] << std::endl;

    w[i] *= 1/integral0;

    integral += w[i];
  }

    
  //  std::cout << "width: " << sigma/n << "\tintegral: " << integral << "\t(" << n << ")" << std::endl;

  //  std::cout << w << "\n" << std::endl;
  
  return w;
}



double fitscale = 1;

