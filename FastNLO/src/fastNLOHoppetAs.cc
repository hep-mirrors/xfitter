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
//  fastNLOAlphas
//  This class inherits the PDF interface from
//  fastNLOLHAPDF, while the alpha_s evolution
//  is superseeded by the Alphas.h class.
//lhasub
//////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <LHAPDF/LHAPDF.h>
//#include "fastnlotk/fastNLOReader.h"
//#include "fastnlotk/speaker.h"
#include "fastNLOLHAPDF.h"
#include "fastNLOHoppetAs.h"
#include "hoppet_v1.h"

using namespace std;



fastNLOHoppetAs::fastNLOHoppetAs(string name) : fastNLOHoppet(name) {
    //Set some meaningful values
    SetPDGValues();
};
fastNLOHoppetAs::fastNLOHoppetAs(string name, string LHAPDFFile, int PDFSet = 0) :
    fastNLOHoppet(name,LHAPDFFile,PDFSet)
    {
        //Set some meaningful values
        SetPDGValues();
    };

vector<double> fastNLOHoppetAs::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   //
   return fastNLOLHAPDF::GetXFX(xp, muf);
}
