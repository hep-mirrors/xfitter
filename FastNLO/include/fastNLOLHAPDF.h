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

#ifndef FASTNLOLHAPDF
#define FASTNLOLHAPDF

#include "fastNLOReader.h"
//#include <iostream>
//#include <cstdio>
//#include <cstdlib>
//#include "speaker.h"

using namespace std;


class fastNLOLHAPDF : public fastNLOReader {

private:
public:
   fastNLOLHAPDF(string name);
   fastNLOLHAPDF(string name, string LHAPDFfile, int PDFSet = 0);

   // Initializer. Necessary for some alternative evolutions.
   virtual void InitEvolveAlphas();
   // Pseudo-Setters. Don´t work with LHAPDF, but print warning instead.
   virtual void SetMz(double Mz);
   virtual void SetNFlavor(int nflavor);
   virtual void SetNLoop(int nloop);
   virtual void SetAlphasMz(double AlphasMz , bool ReCalcCrossSection = false);
   virtual void SetQMass(int pdgid, double mq);
   // Setters
   void SetLHAPDFFilename(string filename);
   void SetLHAPDFMember(int set);
   // Getters
   int GetIPDFMember() const;
   int GetNPDFMembers() const;
   int GetNPDFMaxMember() const;
   void PrintPDFInformation() const ;
   virtual double GetQMass(int pdgid) const;
   int GetNLoop() const;
   int GetNFlavor() const;
   double GetAlphasMz(double Q);

protected:
   // inherited functions
   double EvolveAlphas(double Q) const ;
   bool InitPDF();
   vector<double> GetXFX(double xp, double muf) const ;

   // ---- LHAPDF vars ---- //
   string fLHAPDFFilename;
   #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   LHAPDF::PDFSet* PDFSet;
   LHAPDF::PDF* PDF;
   #endif
   int fnPDFs;
   int fiPDFMember;

   double fchksum;


};

#endif
