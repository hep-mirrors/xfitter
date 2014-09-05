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


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <LHAPDF/LHAPDF.h>
//#include "fastnlotk/speaker.h"
//#include "fastnlotk/fastNLOReader.h"
#include "fastNLOLHAPDF.h"
#include "fastNLOAlphas.h"
#include "Alphas.h"

using namespace std;

fastNLOAlphas::fastNLOAlphas(string name) : fastNLOLHAPDF(name) {
      ;
   };
fastNLOAlphas::fastNLOAlphas(string name, string LHAPDFFile, int PDFSet = 0) : fastNLOLHAPDF(name,LHAPDFFile,PDFSet), fAlphasMz(0.1184) {
      ;
   };

//______________________________________________________________________________
double fastNLOAlphas::GetAlphasMz() const {
      return fAlphasMz;
   };

void fastNLOAlphas::SetMz(double Mz) {
   Alphas::SetMz(Mz);
}

void fastNLOAlphas::SetNFlavor(int nflavor) {
   if (nflavor == 0) {
      Alphas::SetFlavorMatchingOn(true);
      Alphas::SetNf(6);
      warn["SetNFlavor"]<<"GRV evolution of alpha_s is implemented for Nf=5 only.\n";
      warn["SetNFlavor"]<<"You chose a variable Nf with Nfmax=6, i.e. results for Nf other than 5 presumably are wrong!\n";
   } else if (nflavor == 5) {
      Alphas::SetNf(nflavor);
   } else {
      error["SetNFlavor"]<<"GRV evolution of alpha_s is implemented for Nf=5 only.\n";
      exit(1);
   }
}

void fastNLOAlphas::SetNLoop(int nloop) {
   Alphas::SetNLoop(nloop);
}

void fastNLOAlphas::SetAlphasMz(double AlphasMz , bool ReCalcCrossSection) {
   debug["SetAlphasMz"]<<"Setting alpha_s(Mz)="<<AlphasMz<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set the alpha_s value at M_Z
   //
   fAlphasMz    = AlphasMz;             // new alpha_s value
   if (ReCalcCrossSection) CalcCrossSection();
}


//______________________________________________________________________________


double fastNLOAlphas::EvolveAlphas(double Q) const {
   //
   // Implementation of Alpha_s evolution as function of Mu_r only.
   //
   // the alpha_s evolution is done within LHAPDF.
   //
   // WARNING: You cannot change alpha_s(Mz), but is is
   // defined with the pdf. 'alphasMz' is not used here!
   //
   return Alphas::CalcAlphasMu(Q , fAlphasMz);
}


//______________________________________________________________________________


void fastNLOAlphas::SetGRVtoPDG2012_2loop() {
   info["SetGrVtoPDF2012"]<<"Resetting to GRV Alphas::Alphas evolution."<<endl;
   Alphas::SetMz(91.1876); // PDG 2012
   Alphas::SetNf(5);
   Alphas::SetNLoop(2);
   Alphas::SetFlavorMatchingOn(false);
   if (info.GetSpeak()) {
      info<<"Calling Alphas::PrintInfo()."<<endl;
      info<<"Alpha_s(Mz) value is taken from fastNLOAlphas, instead of Alphas::Alphas."<<endl;
      Alphas::PrintInfo();
   }
}
