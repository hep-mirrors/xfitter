#ifndef CIJET_H
#define CIJET_H

//--------------------------------------------------------------
// C++ wraper for using CIJET interface 


#include <vector>
#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include "ReactionTheory.h"
//--------------------------------------------

using namespace std;

class CIJETReader {

 public:

// initialization
  CIJETReader(string fname) {
	filepds=fname;
  }

  CIJETReader(string fname, ReactionTheory* reaction) {
	filepds=fname;
        _reactionTheory=reaction;
  }

// name of table file  
  string filepds;

// QCD parameters
  int Oqcd=1;
  double mur=1.0, muf=1.0;

// contact interaction parameters
  int cshape=1;
  double lambda, c1, c2, c3;

// normalization
  double fut=1.; 

// setting parameters
  void setnorm(double norm);
  void setorder(int od);
  void setscales(double mr, double mf);
  void setci(double lam, vector<double> coe);

// initialization
  void setinit();

// calculating xsecs 
  vector<double> calcxsec();

 private:

// commons
  ReactionTheory* _reactionTheory;

// Get strong coupling
  double GetAs_(double mr);

// Get PDFs (recalculate, e.g. call to xFitter routine)
  void GetPDF_(double x, double mf, double pdf[13]);

};

#endif
