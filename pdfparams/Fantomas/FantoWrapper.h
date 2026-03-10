/*
Description: an extern C wrapper for core Fantomas functions
  suitable for linking to Fortran and other codes
*/

#pragma once

#ifndef FANTOWRAPPER_H
#define FANTOWRAPPER_H

void readfantosteer();
void writefantoout();
void updatefantopars(int &ifl,double *parsin);
double fantopara(int &ifl, double &x);
double fantoMellinMoment(int &ifl, int &MellinPower, int npts=10000);
void getfantochi2(double& fantochi2);

extern "C"
{
  void readfantosteer_();
  
  void writefantoout_();
  
  void updatefantopars_(int &ifl,double *parsin);
  
  double fantopara_(int &ifl, double &x);
  
  double fantomellinmoment_(int &ifl, int &mellinpower, int npts=10000);
  
  void getfantochi2_(double& fantochi2); 

}

#endif //FANTOWRAPPER_H
