// Author: Daniel Britzger
// DESY, 08/08/2012

#ifndef fASTNLODIFFUSER
#define fASTNLODIFFUSER


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  FastNLODiffUSER                                                     //
//                                                                      //
//  FastNLODiffReader is a standalone code for reading                  //
//  diffractive FastNLO tables of version 2.0 for DIS processes         //
//                                                                      //
//  This is a dummy class and not working!!!                            //
//                                                                      //
//  Please insert your desired code into the functions                  //
//   - EvolveAlphas()                                                   //
//   - InitPDF()                                                        //
//   - GetDiffXFX()                                                     //
//                                                                      //
//  Within most applications it is necessary to define the functions    //
//  within a .cc-file.                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <string>
#include <cstdio>
#include <vector>
#include "fastnlotk/fastNLODiffReader.h"


class fastNLODiffUser : public fastNLODiffReader {

public:

   fastNLODiffUser(string filename);
   ~fastNLODiffUser(void) {;};

protected:

   // inherited functions
   double EvolveAlphas(double Q) const ;
   bool InitPDF();
   vector<double> GetDiffXFX(double xpom, double zpom, double muf) const ;

};


//______________________________________________________________________________


#endif
