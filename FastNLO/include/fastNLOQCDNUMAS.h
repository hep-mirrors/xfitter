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
//
//////////////////////////////////////////////////////////////////////////
#ifndef FASTNLOQCDNUMAS
#define FASTNLOQCDNUMAS

//#include "fastNLOReader.h"
//#include <iostream>
//#include <cstdio>
//#include <cstdlib>
//#include <string.h>
//#include <LHAPDF/LHAPDF.h>
//#include "speaker.h"
#include "fastNLOLHAPDF.h"


using namespace std;

extern "C" {
   double asfunc_(double* r2, int* nf  , int* ierr);
   double qcinit_(int* lun, char* filename, int);
   double setalf_(double* alfs, double* r2);
   double setord_(int* iord);
   double setcbt_(int* nfix, int* iqc, int* iqb, int* iqt);
   double gqmake_(double* qarr, double* wgt, int* n, int* nqin, int* nqout);
   int iqfrmq_(double* q2);
}

class fastNLOQCDNUMAS : public fastNLOLHAPDF {

public:
   fastNLOQCDNUMAS(string name);
   fastNLOQCDNUMAS(string name, string LHAPDFFile, int PDFSet);
   //inherited
   void CalcCrossSection();

   void InitEvolveAlphas();
   // ---- Alphas vars ---- //
   // Setters
   void SetMz(double Mz);
   void SetNFlavor(int nflavor);
   void SetNLoop(int nloop);
   void SetQMass(int pdgid, double qmass);
   void SetAlphasMz(double AlphasMz , bool ReCalcCrossSection = false);
   void SetPDGValues();
   void SetLHAPDFValues();
   // Getters
   double GetMz() const;
   double GetQMass(int pdgid) const;
   int GetNFlavor(int nflavor) const;
   int GetNLoop() const;
   double GetAlphasMz() const;



protected:

   // inherited functions
   double EvolveAlphas(double Q) const ;
   // ---- Alphas vars ---- //
   double fAlphasMz;
   double fMz;
   int fnFlavor;
   int fnLoop;
   double QMass[6];


};

#endif
