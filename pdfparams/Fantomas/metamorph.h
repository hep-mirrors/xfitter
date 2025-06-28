/*****************************************************************
 *
 *     Program: Fantomas
 *
 *     Advanced parametrizations of nonpertubative QCD functions
 *
 *     Version 1.0 
 *
 *    Description: the metamorph parametrization class
 *
 * History: June 2025 version 1.0 
 ******************************************************************/
#ifndef METAMORPH_H
#define METAMORPH_H

using namespace std;

#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
// #include <omp.h>
#include <cstddef>
#include <vector>
#include <string>
#include "cl2DArray.h"
#include "LUPinverse.h"

//========================================================================
// Auxiliary functions
//========================================================================

// The default carrier function. The input parameters are the momentum fraction x 
// and a pointer to the array a with control parameters.
double DefCarrier(double x, double *a);

//========================================================================
// Class definitions
//========================================================================

class ControlPoint
//Defines a control point at position x, with the value f returned dependent on 
//the control parameter s in the interval fm <= 0 <= fp. 
{
public:
  double x;  // x position of the control point
  double *ps; // Pointer to an external parameter that will control the value of the function f
  double fm = 0.0, fp = INFINITY; // minimal and maximal values of function f
  //MappingMode determines how the value of f is determined, given s.
  //MappingMode=0: no scaling
  //             1: bounded linear scaling (not implemented)
  //             2: bounded softsign scaling (not implemented)
  int  MappingMode = 0; 

  int SetBoundary(const int MappingModeIn = 0,
                  const double fmIn = 0.0, const double fpIn = INFINITY);
  ControlPoint(const double xIn, double &psIn, const int MappingModeIn = 0,
               const double fmIn = 0.0, const double fpIn = INFINITY);
  double s();
  // the value of the function at the control point, given the parameter s
  double f();
}; //class ControlPoint->

//====================================================
class metamorph
/* metamorph provides a PDF parametrization dependent on the momentum fraction x
 *  of the form
 *       f(x) = Carrier(x, Sc)*Modulator(x, Sm),
 *  where 
 *   Carrier defines the asymptotic behavior of f(x) in the limits x->0 and x->1;
 *   Modulator(x,Sm) determines the behavior of f(x) over the interval 0 < x < 1.
 *
 *    Carrier(x, Sc) is a smooth function dependent on the array of external 
 *     parameters Sc[], such as Carrier = x^Sc(0)*(1-x)^Sc(1).
 *   Modulator is given by a Bezier curve of order Ncp-1 that depends on x^xPower 
 *     as its argument, where default xPower=0.5.
 *     The shape of Modulator(x, Sm) is controlled by external parameters Sm[] 
 *     that set the value of f(x) at Ncp user-selected control points CP at discrete x values 
 *     Xs[0], Xs[1], ..., Xs[Ncp-1]. Internally, Modulator(x, Sm) is controlled by a locally stored
 *     vector of Bezier coefficients C. 
 *   UpdateModulator() must be called every time the external parameters Sm are 
 *     changed to update the vector C before the function f(x) is called. 
 */
{
  int Nm; // numbers of free parameters in the carrier and modulator functions
  vector<ControlPoint> CP; // Control points for the modulator function
  cl2DArray<double> *M = NULL, *InvM = NULL; // Matrices M and T for the modulator function, and their inverse matrices
  cl2DArray<double> *T = NULL, *InvT = NULL, *multi = NULL;
  double *C, *P;     // Bezier coefficients of the modulator function and a vector of corresponding parameters
  int *Pdec;
  
  double (*Carrier)(double, double*); // a pointer to the carrier function. Points to function DefCarrier by default. 
  bool ModulatorSet = false;

  double  xPower=0.5, xstrmin = 1e-10, xstrmax = 1.0-1e-10,
    stretchPower=6.0;        //parameters of the x stretching function y(x)
                             //in the Bezier curve, their default values

public:
  double *pSc;    // pointer to an external array with parameters of the modulator function
  int ID=999;     //the ID of the metamorph in the collection, set externally  
  
  metamorph(const int NmIn, const double XsIn[], double SmIn[], double ScIn[], const double xPowerIn = 0.5,
            const vector<double>& vstretch={}, double (*CarrierIn)(double, double *) = &DefCarrier);
  // Constructor: construct a metamorph class by reading Nm values from an 
  // input array Xs and setting control point parameters to point to an array 
  // of external variables Sm. Also provide a pointer to the array Sc containing the parameters
  // for the carrier function. DefCarrier is the default carrier function.
  // An optional input parameter CarrierIn is a pointer to an alternative carrier function.
  // An optional vector vtretch passes xstrmin and xstrmax values for the stretching function

  int SetBoundary(const int MappingModeIn, double fmIn[], double fpIn[]);
  //tbd Set the boundaries of the possible modulator values

  void UpdateModulator();                    
  // Compute matrix C, given the carrier function vectors Fc, matrices InvM and InvT. This computation 
  // is done every time the control points are updated, and before the metamorph value f is computed
  // for varied x values. We write the PDF function as f(x_i) = Fc(x_i) * Fm(x_i),
  // where Fc is the carrier, Fm is the modulator.
  // A control point returns f_i with f_-i <= f_i <= f_+i. Then we have f_i = Fc_i * Fm_i. 
  // The vector of Bezier coefficients C is computed as C = InvM * InvT * P, where P 
  // is the vector of P_i = Fm(x_i).

  double f(const double x);
  // Returns the value of the metamorph for the momentum fraction x

  //tbd void SetModulatorXPower(double xPowerIn); //Set xPower, recompute T if necessary
  //tbd  SetCustomCarrier(double *CarrierIn); //Set a pointer to a custom carrier function
  //tbd double GetModulatorXPower();    //Get the value of xPower

  double GetMellinMoment(double MellinPower, int npts = 10000);
  // metamorph::GetMellinMoment() returns <x^(n+1) f>, i.e., the  integral of x^(n+1) * f(x) * dx
  // over 0 < x < 1, where n=MellinPower.

  double GetConditionNumber();
  // Returns the condition number of the matrix T and T^{-1} using the 
  // Frobenius norm, ||T||, i.e. the square root of the sum of absolute
  //squares of the elements. The condition number is calculated
  //as ||T|| * ||T^{-1}||. 

  double Cs(const int i);
  // Returns the value for the ith Bezier coefficient

  double Modulator(const double x);
  // Returns the value of the Modulator function for the momentum fraction x

  double yx(const double x);
  //Returns the stretched argument y of the Bezier curve, y(x)=xstretched(x)^Npower 
  
  void SetXstretching(const double xstrminIn=1e-10, const double xstrmaxIn=1.0-1e-10, const double stretchPowerIn=6.0);
  //Sets parameters of the stretching function yx

  void GetXstretching(double xstrminOut, double xstrmaxOut, double stretchPowerOut);
  //Sets parameters of the stretching function yx
 
 
  ~metamorph();
}; //class metamorph -> 

#endif //METAMORPH_H
