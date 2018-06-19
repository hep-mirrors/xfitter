// emacs: this is -*- c++ -*-
//
//   appl_pdf.h        
//
//   pdf transform functions header                  
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: appl_pdf.h, v   Fri Dec 21 22:19:50 GMT 2007 sutt $


#ifndef __APPL_PDF_H
#define __APPL_PDF_H

#include <iostream>
#include <fstream>
#include <sstream>

#include <vector> 
#include <map> 
#include <string> 

#include <exception> 


namespace appl { 


class appl_pdf;

typedef std::map<const std::string, appl_pdf*> pdfmap;


// this is a *maybe* nice class, a base class for pdf 
// functions
//
// it has a virtual evaluate() method to be definied in 
// the derived class, and a static std::map of all the names
// of instances of the derived classes
//
// when a new instance of the class is created, it 
// automatically adds it's name to the std::map, so the user 
// doesn't need to worry about consistency, and removes 
// itself when the derived instance is deleted

class appl_pdf { 

public:

  // pdf error exception
  class exception : public std::exception { 
  public: 
    exception(const std::string& s="") { std::cerr << what() << " " << s << std::endl; }; 
    //exception(std::ostream& s)         { std::cerr << what() << " " << s << std::endl; }; 
    exception(std::ostream& s)         { std::stringstream ss; ss << s.rdbuf(); std::cerr << what() << " " << ss.str() << std::endl; }; 
    const char* what() const throw() { return "appl::appl_pdf::exception "; }
  };
  
public:

  /// constructor and destructor
  appl_pdf(const std::string& name);

  virtual ~appl_pdf();

  /// retrieve an instance from the std::map 
  static appl_pdf* getpdf(const std::string& s, bool printout=true);
  
  /// print out the pdf std::map
  static void printmap(std::ostream& s=std::cout) {
    pdfmap::iterator itr = __pdfmap.begin();
    while ( itr!=__pdfmap.end() )  {
      s << "pdfmap " << itr->first << "\t\t" << itr->second << std::endl;
      itr++;
    } 
  }

  /// initialise the factory  
  static bool create_map(); 

  virtual void evaluate(const double* fA, const double* fB, double* H) = 0; 

  virtual int decideSubProcess( const int , const int  ) const;

  std::string   name() const { return m_name;  }

  int     Nproc() const { return m_Nproc; } 
  int     size()  const { return m_Nproc; } 



  std::string  rename(const std::string& name) { 
    /// remove my entry from the std::map, and add me again with my new name
    if ( __pdfmap.find(m_name)!=__pdfmap.end() ) { 
      __pdfmap.erase(__pdfmap.find(m_name));
    }
    else { 
      std::cout << "appl_pdf::rename() " << m_name << " not in std::map" << std::endl;
    }
    m_name = name;
    addtopdfmap(m_name, this);
    return m_name;
  }


  /// code to allow optional std::vector of subprocess contribution names

  const std::vector<std::string>& subnames() const { return m_subnames; }

  void addSubnames( const std::vector<std::string>& subnames ) { m_subnames = subnames; }

  void  addSubname( const std::string& subname ) { 
    if ( int(m_subnames.size())<m_Nproc-1 ) m_subnames.push_back(subname); 
  }



  /// is this a W+ or a W- pdf combination? or neither?
  int getckmcharge() const { return m_ckmcharge; }

  /// access the ckm matrices - if no matrices are required these std::vectors have 
  /// zero size

  const std::vector<double>&               getckmsum() const { return m_ckmsum; }
  const std::vector<std::vector<double> >& getckm2()   const { return m_ckm2; }
  const std::vector<std::vector<double> >& getckm()    const { return m_ckm; }
  
  /// set the ckm matrices from external values
  void setckm( const std::vector<std::vector<double> >& ckm ); 
  void setckm2( const std::vector<std::vector<double> >& ckm2 ); 

  /// code to create the ckm matrices using the hardcoded default 
  /// values if required
  /// takes a bool input - if true creates the ckm for Wplus, 
  /// false for Wminus
  void make_ckm( bool Wp=true );

  void SetNProc(int Nsub){ m_Nproc=Nsub; return;};

  /// set some useful names for the different subprocesses
  void setnames( const std::vector<std::string>& names) { m_names = names; } 
  std::vector<std::string> getnames() const  { return m_names; } 
  
protected:

  /// search the path for configuration files
  static  std::ifstream& openpdf( const std::string& filename ); 

private:

  static void addtopdfmap(const std::string& s, appl_pdf* f) { 
    if ( __pdfmap.find(s)==__pdfmap.end() ) { 
      __pdfmap.insert( pdfmap::value_type( s, f ) );
      //      std::cout << "appl_pdf::addtomap() registering " << s << " in std::map addr \t" << f << std::endl;
    }
    else { 
      throw exception( std::cerr << "appl_pdf::addtopdfmap() " << s << " already in std::map\t0x" << __pdfmap.find(s)->second  );
    }
  }
  
protected:

  int         m_Nproc;
  std::string m_name;

  std::vector<std::string> m_subnames;

  /// ckm matrix related information 
  /// W+, W- or neither?
  int  m_ckmcharge;

  // ckm matrices
  std::vector<double>               m_ckmsum;
  std::vector<std::vector<double> > m_ckm2; /// squared 13x13 matrix
  std::vector<std::vector<double> > m_ckm;  /// simple 3x3

  /// some strings for more useful name if required
  std::vector<std::string>           m_names;

  static pdfmap                     __pdfmap;
  static std::vector<std::string>   __pdfpath;
};


}


#endif  // __APPL_PDF_H 










