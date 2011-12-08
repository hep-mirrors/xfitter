#include <math.h>
#include "xsect_dipole_proton.h"

/////////////////////////////////////////////////////////
// Parametrisation of the dipole-proton cross-section  //
// as obtained from high-energy QCD with saturation    //
/////////////////////////////////////////////////////////
// This file contains two version of the cross-section //
// The first one is the original from 		       //
//   [1] E. Iancu, K. Itakura and S. Munier, 	       //
//       Phys.Lett.B590:199-208,2004 [hep-ph/0310338]  //
// which includes only light quarks.		       //
// The second is from				       //
//   [2] G. Soyez, in preparation		       //
// which also accounts for heavy quarks.	       //
// They both implement eq. (8) of Ref. 1 or eq. (3) of //
// Ref. 2. Only the parameters of the fit vary.	       //
//						       //
// - Arguments:					       //
//    r : dipole size				       //
//    Y : rapidity (Y=log(1/x))			       //
// - Returned value:				       //
//    the dipole-proton cross-section (in GeV^{-2})    //
// 						       //
// Additional remarks:				       //
// - The quark masses are			       //
//     m_{u,d,s}=140 MeV, m_c=1.4 GeV and m_b=4.5 GeV  //
// - When the heavy quarks are considered (2nd model), //
//   the x value for their contribution has to be      //
//   shifted to x(1+4 m^2/Q^2)			       //
/////////////////////////////////////////////////////////

// light quarks only
//-------------------
double sigma_iim_(double *rr, double *YY
		  , double *lamf, double *x0f, double *rf){

  double r = *rr;
  double Y = *YY;

  if (Y<=0) return 0.0;

  // parameters (common to the 2 models)
  double N0 = 0.7;
  double kappa = 9.9;

  // parameters (specific to the light-flavour fit
  double gc = 0.6275;

  // These are fit parameters:

  double lambda = *lamf; // 0.253;
  double x0 = *x0f; // 0.267e-4;
  double Rp2 =*rf * (*rf); // 3.25*3.25;

  // saturation scale and geometric scaling variable
  double Qs2 = exp(lambda*(Y+log(x0)));
  double lQs = 0.5*log(Qs2);
  double z = -log(r/2.0)-lQs;
  
  // McLerran-Venugopalan parameters
  double a = gc*N0/(1-N0);
  a = a*a/log(1/(1-N0));
  double lrqsb = z-(1-N0)/gc/N0*log(1/(1-N0));

  // amplitude
  return 2.0*M_PI*Rp2*
    ((z>0) ? N0*exp(-2.0*z*(gc+z/(kappa*lambda*Y))) 
           : 1-exp(-a*lrqsb*lrqsb));
}

// light+heavy quarks
//--------------------
double sigma_proton_(double *rr, double *YY){

  double r = *rr;
  double Y = *YY;

  if (Y<=0) return 0.0;

  // parameters (common to the 2 models)
  double N0 = 0.7;
  double kappa = 9.9;

  // parameters (specific to the light-flavour fit
  double gc = 0.7376;
  double lambda = 0.2197;
  double x0 = 0.1632e-4;
  double Rp2 = 3.344*3.344;

  // saturation scale and geometric scaling variable
  double Qs2 = exp(lambda*(Y+log(x0)));
  double lQs = 0.5*log(Qs2);
  double z = -log(r/2.0)-lQs;
  
  // McLerran-Venugopalan parameters
  double a = gc*N0/(1-N0);
  a = a*a/log(1/(1-N0));
  double lrqsb = z-(1-N0)/gc/N0*log(1/(1-N0));

  // amplitude
  return 2.0*M_PI*Rp2*
    ((z>0) ? N0*exp(-2.0*z*(gc+z/(kappa*lambda*Y))) 
           : 1-exp(-a*lrqsb*lrqsb));
}

