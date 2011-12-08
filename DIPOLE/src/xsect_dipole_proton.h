#ifndef __XSECT_DIPOLE_PROTON_H__
#define __XSECT_DIPOLE_PROTON_H__

// 25/07/2010
// Added "_" to the names of the procedures, so
// the names have two "_" after last letter.

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
double sigma_iim__(double *r, double *Y
 , double *lamf, double *x0f, double *rf);

// light+heavy quarks
//--------------------
double sigma_proton__(double *r, double *Y);

#endif
