#ifndef FANTOMAS_H
#define FANTOMAS_H

void readfantosteer();
void writefantoout();
void updatefantopars(int &flavor,double *parsin);
double fantopara(int &flavor, double &x);
double fantoMellinMoment(int &flavor, int &MellinPower, int npts=10000);
void getfantochi2(double& fantochi2);
#endif
