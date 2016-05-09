#ifndef CUTE_H
#define CUTE_H

//Cute header
#include "plot.h"

#include <string>
#include <vector>


using namespace std;

class Cute 
{
 public:
  Cute(string file);

  void SetBins(vector <double> lowedge, vector <double> upedge, int bindensity);

  void Calculate(const double mur, const double muf);
  
  //values of the cross section in each bin
  vector<double> values;
  //bin widths
  vector<double> binwidth;

  //points where the cross section is evaluated, map xs into qt
  map <double, double> points;

  //points-per-bin: steer this parameter
  int bindensity;

  string infile;

  //nominal accuracy, as set in the input file
  double accuracy;

  Plot plot;
};


#endif
