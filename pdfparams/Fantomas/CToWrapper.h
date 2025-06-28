/*
Description: an extern C wrapper for core Fantomas functions
  suitable for linking to Fortran and other codes
*/

#ifndef CTOWRAPPER_H
#define CTOWRAPPER_H

extern "C"{
  void readfantosteer_();
  void writefantoout_();
  void updatefantopars_(int &flavor,double *parsin);
  double fantopara_(int &flavor, double &x);
  double fantomellinmoment_(int &flavor, int &MellinPower, int npts=10000);
  //void getfantochi2(double& fantochi2);
}
#endif //CTOWRAPPER_H
