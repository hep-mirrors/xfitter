// Author: Daniel Britzger
// DESY, 20/04/2012

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  fastNLO_reader_2.1.0                                                //
//  D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch           //
//                                                                      //
//  The projects web page can be found at:                              //
//    http://projects.hepforge.org/fastnlo                              //
//                                                                      //
//  If you use this code, please cite:                                  //
//    T. Kluge, K. Rabbertz and M. Wobisch, hep-ph/0609285              //
//    D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch,        //
//       arXiv:1109.1310                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//
//  fastNLOAlhpas
//  This class inherits the PDF interface from
//  fastNLOLHAPDF, while the alpha_s evolution
//  is superseeded by the Alphas.h class.
//
//////////////////////////////////////////////////////////////////////////


#ifndef FASTNLOALPHAS
#define FASTNLOALPHAS

//#include "fastNLOReader.h"
//#include <iostream>
//#include <cstdio>
//#include <cstdlib>
//#include <LHAPDF/LHAPDF.h>
//#include "speaker.h"
#include "fastNLOLHAPDF.h"
//#include "Alphas.h"

using namespace std;

class fastNLOAlphas : public fastNLOLHAPDF {

public:
   fastNLOAlphas(string name);
   fastNLOAlphas(string name, string LHAPDFFile, int PDFSet);

   // ---- Alphas vars ---- //
   // Setters
   void SetMz(double Mz);
   void SetNFlavor(int nflavor);
   void SetNLoop(int nloop);
   void SetAlphasMz(double AlphasMz , bool ReCalcCrossSection = false);
   // Getters
   double GetAlphasMz() const;
   void SetGRVtoPDG2012_2loop();


protected:
   // inherited functions
   double EvolveAlphas(double Q) const ;

   // ---- Alphas vars ---- //
   double fAlphasMz;

};



//______________________________________________________________________________


#endif
