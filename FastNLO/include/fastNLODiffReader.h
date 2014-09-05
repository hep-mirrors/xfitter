// Author: Daniel Britzger
// DESY, 02/04/2012

#ifndef fASTNLODIFFREADER
#define fASTNLODIFFREADER


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  fastNLODiffReader                                                   //
//                                                                      //
//  fastNLODiffReader is a standalone code for reading                  //
//  diffractive FastNLO tables of version 2.0 for DIS processes         //
//                                                                      //
//    http://projects.hepforge.org/fastnlo                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <string>
#include <vector>
#include "fastNLOReader.h"

using namespace std;

class fastNLODiffReader : public fastNLOReader {

public:

   fastNLODiffReader(string filename);
   virtual ~fastNLODiffReader(void) {};

   void SetXPomSlicing(int nStep, double* xpom, double* dxpom);
   void SetXPomLogSlicing(int nStep, double xpommin, double xpommax);
   void SetXPomLinSlicing(int nStep, double xpommin, double xpommax);
   void SetXPomExpSlicing(int nStep, double xpommin, double xpommax);
   void SetZRange(double zmin , double zmax) {
      fzmin = zmin ;
      fzmax = zmax;
   };
   double GetZRangeMin() {
      return fzmin;
   };
   double GetZRangeMax() {
      return fzmax;
   };

   vector < double > GetCrossSection();
   void CalcCrossSection();
   vector < double > GetDiffCrossSection();
   void FillPDFCache(bool ReCalcCrossSection = false);
   vector < double > GetReferenceCrossSection();

   // ---- Print outs must be overwritten ---- //
   void PrintCrossSectionsWithReference();


protected:

   double fxpom;
   double fzmin;
   double fzmax;

   vector < double > fxPoms;
   vector < double > fdxPoms;

   // inherited functions
   virtual double EvolveAlphas(double Q) const = 0;
   virtual bool InitPDF() = 0;
   vector<double> GetXFX(double xp, double muf) const;
   virtual vector<double> GetDiffXFX(double xpom, double zpom, double muf) const = 0;

};

#endif

