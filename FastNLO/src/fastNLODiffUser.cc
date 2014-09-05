// Author: Daniel Britzger
// DESY, 08/08/2012


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  fastNLODiffUSER                                                     //
//                                                                      //
//  fastNLODiffReader is a standalone code for reading                  //
//  diffractive fastNLO tables of version 2.0 for DIS processes         //
//                                                                      //
//  This is a dummy class and not working!!!                            //
//                                                                      //
//  Please insert your desired code into the functions                  //
//   - EvolveAlphas()                                                   //
//   - InitPDF()                                                        //
//   - GetDiffXFX()                                                     //
//                                                                      //
//  Within most applications it is necessary to define the functions    //
//  within a .cc-file explicitly.                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <string>
#include <cstdio>
#include <vector>
#include "fastNLODiffUser.h"


fastNLODiffUser::fastNLODiffUser(string filename) : fastNLODiffReader(filename) {
}


//______________________________________________________________________________


double fastNLODiffUser::EvolveAlphas(double Q) const {
   // --- fastNLO user:
   // Implementation of Alpha_s evolution as function of the
   // factorization scale [and alphas(Mz)].
   //
   return 0;
}


//______________________________________________________________________________


bool fastNLODiffUser::InitPDF() {
   // --- fastNLO user:
   //  Initalize PDF parameters if necessary
   //
   return true;
}


//______________________________________________________________________________



vector<double> fastNLODiffUser::GetDiffXFX(double xpom, double zpom, double muf) const {
   //
   //  GetDiffXFX is used to get the parton array from the
   //  pdf-interface. It should return a vector of 13
   //  parton flavors from tbar to t at a certain
   //  xpom, zpom and factorisation scale.
   //
   vector < double > xfx(13);
   // fastNLO user:
   //   include some function here to fill the parton density array
   //   xfx[0]=tbar, xfx[6]=gluon, xfx[12]=t
   //debug<<"xpom="<<xpom<<"\tzpom="<<zpom<<"\tmuf="<<muf<<"\tgluon = "<<xfx[6]<<endl;
   return xfx;
}


//______________________________________________________________________________

