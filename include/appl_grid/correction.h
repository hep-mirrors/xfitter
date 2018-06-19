// emacs: this is -*- c++ -*-
//
//   @file    correction.h        
//            class to store the multipliciative post processing
//            corrections to be applied, only basic for the time 
//            but will be extended as appropriate
//  
//   Copyright (C) 2014 M.Sutton (sutt@cern.ch)    
//
//   $Id: correction.h, v0.0   Sun 23 Mar 2014 09:08:46 CET sutt $


#ifndef  CORRECTION_H
#define  CORRECTION_H

#include <iostream>
#include <vector>
#include <string>


// typedef std::vector<double> correction;


class correction {

public:

  correction(const std::vector<double>& v, const std::string& s="" ) : mlabel(s), mv(v) { } 

  virtual ~correction() { } 

  std::string label() const { return mlabel; }

  unsigned size() const { return mv.size(); }

  double& operator[](int i)       { return mv[i]; }
  double  operator[](int i) const { return mv[i]; }

  operator std::vector<double>&() { return mv; } 

  correction operator=(const std::vector<double>& v) { mv=v; return *this; } 

private:

  std::string         mlabel;   
  std::vector<double> mv;  

};


// inline std::ostream& operator<<( std::ostream& s, const correction& /* _c */ ) { 
//   return s;
// }



#endif  // CORRECTION_H 










