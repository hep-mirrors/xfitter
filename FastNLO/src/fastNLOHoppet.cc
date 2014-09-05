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
#include "fastNLOHoppet.h"
#include "HoppetInterface.h"

using namespace std;



fastNLOHoppet::fastNLOHoppet(string name) : fastNLOLHAPDF(name) {
    //Set some meaningful values
    SetLHAPDFValues();
};
fastNLOHoppet::fastNLOHoppet(string name, string LHAPDFFile, int PDFSet = 0) :
    fastNLOLHAPDF(name,LHAPDFFile,PDFSet)
    {
        //Set some meaningful values
        SetLHAPDFValues();
    };
// Getters
double fastNLOHoppet::GetMz() const {
    return HoppetInterface::fMz;
}
double fastNLOHoppet::GetQMass(int pdgid) const {
    return HoppetInterface::QMass[pdgid];
}
int fastNLOHoppet::GetNFlavor() const {
    return HoppetInterface::fnFlavor;
}
int fastNLOHoppet::GetNLoop() const {
    return HoppetInterface::fnLoop;
}
double fastNLOHoppet::GetAlphasMz() const {
    return HoppetInterface::fAlphasMz;
};


void fastNLOHoppet::SetPDGValues() {
   // Initialize with PDG values
   HoppetInterface::QMass[0]  = PDG_MD;
   HoppetInterface::QMass[1]  = PDG_MU;
   HoppetInterface::QMass[2]  = PDG_MS;
   HoppetInterface::QMass[3]  = PDG_MC;
   HoppetInterface::QMass[4]  = PDG_MB;
   HoppetInterface::QMass[5]  = PDG_MT;
   HoppetInterface::fMz       = PDG_MZ;
   HoppetInterface::fAlphasMz = PDG_ASMZ;
   //Variable flavor number scheme
   HoppetInterface::fnFlavor = 5;
   //2-loop alpha_s evolution
   HoppetInterface::fnLoop = 2;
   HoppetInterface::InitHoppet(*this);
}

void fastNLOHoppet::SetLHAPDFValues() {
   //Be sure LHAPDF is initialized when reading the properties
   if (fchksum == 0 || fchksum != CalcChecksum(1.)) {
      InitPDF();
   }
   //How to read LHAPDF Mz???
   HoppetInterface::fMz = PDG_MZ;
   HoppetInterface::fAlphasMz = LHAPDF::alphasPDF(HoppetInterface::fMz);
   HoppetInterface::fnLoop = LHAPDF::getOrderAlphaS();
   HoppetInterface::fnFlavor = LHAPDF::getNf();
   for (int i = 0; i < 6; i++)
      HoppetInterface::QMass[i] = LHAPDF::getQMass(i+1);
   HoppetInterface::InitHoppet(*this);
}

void fastNLOHoppet::SetMz(double Mz) {
   HoppetInterface::fMz = Mz;
   HoppetInterface::InitHoppet(*this);
}

void fastNLOHoppet::SetNFlavor(int nflavor) {
   HoppetInterface::fnFlavor = nflavor;
   HoppetInterface::InitHoppet(*this);
}

void fastNLOHoppet::SetNLoop(int  nloop) {
   HoppetInterface::fnLoop = nloop;
   HoppetInterface::InitHoppet(*this);
}

void fastNLOHoppet::SetQMass(int pdgid, double qmass) {
   HoppetInterface::QMass[pdgid] = qmass;
   HoppetInterface::InitHoppet(*this);
}

void fastNLOHoppet::SetAlphasMz(double AlphasMz) {

   HoppetInterface::fAlphasMz    = AlphasMz;             // new alpha_s value
   HoppetInterface::InitHoppet(*this);
}

double fastNLOHoppet::EvolveAlphas(double Q ) const {
   //return HoppetInterface::EvolveAlphas(Q);
   return LHAPDF::alphasPDF(Q);
}

vector<double> fastNLOHoppet::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   //
   return HoppetInterface::GetXFX(xp, muf);
}
