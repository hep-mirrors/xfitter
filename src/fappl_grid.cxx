//   
//   fappl_grid.cxx        
//     
//      fortran callable wrapper functions for the c++  
//      appl grid project.
//                   
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: fappl.cxx, v1.0   Wed May 21 14:31:36 CEST 2008 sutt
//   
//   08/04/2011 AS: changed names to include in xFitter

#include <map>
#include <iostream>
#include <sstream>
#include <string>
using std::cout;
using std::endl;
using namespace std;

#include "xfitter_cpp.h"

#include "appl_grid/appl_grid.h"
#include "appl_grid/fastnlo.h"

/// create a grid
extern "C" void appl_bookgrid_(int& id, const int& Nobs, const double* binlims);

/// delete a grid
extern "C" void appl_releasegrid_(int& id);

/// delete all grids
extern "C" void appl_releasegrids_();

/// read a grid from a file
extern "C" void appl_readgrid_(int& id, const char* s);

/// write to a file 
extern "C" void appl_writegrid_(int& id, const char* s);

/// add an entry 
extern "C" void appl_fillgrid_(int& id, 
			  const int& ix1, const int& ix2, const int& iQ, 
			  const int& iobs, 
			  const double* w,
			  const int& iorder );  

/// redefine the grid dimensions
extern "C" void appl_redefine_(int& id, 
			  const int& iobs, const int& iorder, 
			  const int& NQ2, const double& Q2min, const double& Q2max, 
			  const int& Nx,  const double&  xmin, const double&  xmax); 

/// get number of observable bins for a grid 
extern "C" int appl_getnbins_(int& id);


/// get bin number for a grid
extern "C" int appl_getbinnumber_(int& id, double& data);

/// get low edge for a grid/bin number
extern "C" double appl_getbinlowedge_(int& id, int& bin);

/// get bin width for a grid/bin number
extern "C" double appl_getbinwidth_(int& id, int& bin);

/// do the convolution!! hooray!!
extern "C" void appl_convoluteOld_(int& id, double* data);
extern "C" void appl_convoluteorder_(int& id, int& nloops, double& muR, double& muF, double* data);

extern "C" void appl_convolutewrap_(int& id, double* data, 
			       void (*pdf)(const double& , const double&, double* ),
			       double (*alphas)(const double& ) );


/// print a grid
extern "C" void appl_printgrid_(int& id);

/// print all grids
extern "C" void appl_printgrids_();

/// print the grid documentation
extern "C" void appl_printgriddoc_(int& id);

/// create grids from fastnlo
extern "C" void appl_readfastnlogrids_(  int* ids, const char* s );


static int idcounter = 0;
static std::map<int,appl::grid*> _grid;


/// grid map management

extern "C" void appl_ngrids_(int& n) { n=_grid.size(); }

extern "C" void appl_gridids_(int* ids) { 
  std::map<int,appl::grid*>::iterator gitr = _grid.begin();
  for ( int i=0 ; gitr!=_grid.end() ; gitr++, i++ ) ids[i] = gitr->first;
}



void appl_bookgrid_(int& id, const int& Nobs, const double* binlims) 
{
  id = idcounter++;

  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);

  if ( gitr==_grid.end() ) {
    cout << "bookgrid_() creating grid with id " << id << endl; 
    _grid.insert(  std::map<int,appl::grid*>::value_type( id, new appl::grid( Nobs, binlims,
									      2,    10, 1000, 1,
									      12,  1e-5, 1, 3, 
									      "nlojet", 1, 3, "f3") ) ) ;									 
    //  _grid->symmetrise(true);
  }
  else throw appl::grid::exception( std::cerr << "grid with id " << id << " already exists" << std::endl );  

}


void appl_readgrid_(int& id, const char* s) {
  // fix string from fortran input -- AS
  stringstream istr(stringstream::in | stringstream::out);
  istr << string(s);
  string sstr;
  istr >> sstr;
  id = idcounter++;
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr==_grid.end() ) { 
    _grid.insert(  std::map<int,appl::grid*>::value_type( id, new appl::grid(sstr.c_str()) ) );
  }
  else throw appl::grid::exception( std::cerr << "grid with id " << id << " already exists" << std::endl );  
}


  
void appl_printgrid_(int& id) { 
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) { 
    std::cout << "grid id " << id << "\n" << *gitr->second << std::endl;
  }
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
}


void appl_printgrids_() { 
  std::map<int,appl::grid*>::iterator gitr = _grid.begin();
  for ( ; gitr!=_grid.end() ; gitr++ ) { 
    std::cout << "grid id " << gitr->first << "\n" << *gitr->second << std::endl;
  }
}

  
void appl_printgriddoc_(int& id) { 
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) { 
    std::cout << "grid id " << id << "\n" << gitr->second->getDocumentation() << std::endl;
  }
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
}


void appl_releasegrid_(int& id) { 
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() )     { 
    delete gitr->second; 
    _grid.erase(gitr);
  }
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
}


void appl_releasegrids_() { 
  std::map<int,appl::grid*>::iterator gitr = _grid.begin();
  for ( ; gitr!=_grid.end() ; gitr++ ) { 
    delete gitr->second; 
    _grid.erase(gitr);
  }
}


void appl_redefine_(int& id, 
	       const int& iobs, const int& iorder, 
	       const int& NQ2, const double& Q2min, const double& Q2max, 
	       const int& Nx,  const double&  xmin, const double&  xmax) 
{
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) {
    gitr->second->redefine(iobs, iorder, 
			   NQ2, Q2min, Q2max, 
			   Nx,   xmin,  xmax); 
  }
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
  
} 



int appl_getnbins_(int& id) { 
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) return gitr->second->Nobs();
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
}

int appl_getbinnumber_(int& id, double& data) { 
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) return gitr->second->obsbin(data);
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
}

double appl_getbinlowedge_(int& id, int& bin) { 
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) return gitr->second->obslow(bin);
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
}

double appl_getbinwidth_(int& id, int& bin) { 
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) return gitr->second->deltaobs(bin);
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
}


void appl_convoluteOld_(int& id, double* data) { 
  appl_convolutewrap_(id, data, appl_fnpdf_, appl_fnalphas_); 
}


void appl_convoluteorder_(int& id, int& nloops, double& muR, double& muF, double* data) { 
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) { 
    appl::grid*    g = gitr->second;
    vector<double> v = g->vconvolute(appl_fnpdf_, appl_fnalphas_, nloops, muR, muF);
    for ( unsigned i=0 ; i<v.size() ; i++ ) data[i] = v[i];      
  }
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
}

void appl_convolutewrap_(int& id, double* data, 
		       void (*pdf)(const double& , const double&, double* ),  
		       double (*alphas)(const double& ) ) {  
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) { 
    appl::grid*    g = gitr->second;
    vector<double> v = g->vconvolute( pdf, alphas);
    for ( unsigned i=0 ; i<v.size() ; i++ ) data[i] = v[i];      
  }
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
}




void _appl_writegrid_(int& id, const char* s) { 
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) { 
    cout << "writegrid_() writing " << s << "\tid " << id << endl;
    appl::grid* g = gitr->second;
    g->trim();
    //   g->print();
    g->Write(s);
  }
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
}



void appl_fillgrid_(int& id, 
	       const int& ix1, const int& ix2, const int& iQ,  
	       const int& iobs, 
	       const double* w, 
	       const int& iorder ) { 
  //  cout << "ix " << ix1 << " " << ix2 << "  iQ" << iQ << " " << iobs << "  " << iorder << endl;  
  std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
  if ( gitr!=_grid.end() ) { 
    gitr->second->fill_index(ix1, ix2, iQ, iobs, w, iorder);
  }  
  else throw appl::grid::exception( std::cerr << "No grid with id " << id << std::endl );
 
}


void appl_readfastnlogrids_( int* ids, const char* s ) { 

  /// create the fastnlo grids
  fastnlo f(s);

  /// don't want the grids managed by the fastnlo object, 
  /// manage them in fortran with the map
  f.manageGrids(false);

  ///copy to the fortran accessible grid map
  std::vector<appl::grid*> grids = f.grids();

  //  std::cout << "hooray!" << std::endl;
  
  for ( unsigned i=0 ; i<grids.size() ; i++ ) { 
    int id = idcounter++;
    std::map<int,appl::grid*>::iterator gitr = _grid.find(id);
    if ( gitr==_grid.end() )  { 
      _grid.insert(  std::map<int,appl::grid*>::value_type( id, grids[i] ) );
      // std::cout << grids[i]->getDocumentation() << std::endl;
    }
    else throw appl::grid::exception( std::cerr << "grid with id " << id << " already exists" << std::endl );
    ids[i] = id;
  }  

}

  

