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
#include "fastNLOConstants.h"
#include <LHAPDF/LHAPDF.h>
#include <cmath>


class fastNLOLHAPDF : public fastNLOReader {

private:
public:
   fastNLOLHAPDF(std::string name);
   fastNLOLHAPDF(const fastNLOTable&);
   ~fastNLOLHAPDF();
   fastNLOLHAPDF(std::string name, std::string LHAPDFfile, int PDFSet = 0);
   fastNLOLHAPDF(const fastNLOTable&, std::string LHAPDFfile, int PDFSet = 0);

   // Initializer. Necessary for some alternative evolutions.
   virtual void InitEvolveAlphas();
   // Pseudo-Setters. Don´t work with LHAPDF, but print warning instead.
   virtual void SetMz(double Mz);
   virtual void SetNFlavor(int nflavor);
   virtual void SetNLoop(int nloop);
   virtual void SetAlphasMz(double AlphasMz);
   virtual void SetQMass(int pdgid, double mq);
   // Setters
   void SetLHAPDFFilename(std::string filename);
   void SetLHAPDFMember(int set);
   // Getters
   std::string GetLHAPDFFilename() const {return fLHAPDFFilename;}
   int GetIPDFMember() const;
   int GetNPDFMembers() const;
   int GetNPDFMaxMember() const;
   void PrintPDFInformation() const ;
   virtual double GetQMass(int pdgid) const;
   int GetNLoop() const;
   int GetNFlavor() const;
   LHAPDF::PDFSet* GetPDFSet() const { return PDFSet;};
   LHAPDF::PDF* GetPDF() const { return PDF;};

   double GetAlphasMz() const;

   //! Return struct with vectors (for C++) or vector of vectors (for Python) containing the cross section values and the selected uncertainty
   //! Enum of Uncertaintstyle decides on method to call, but does not work for Python extension --> switch back to use differently named UncertaintyVec methods
   // Use implementations in fastNLOReader for these
   XsUncertainty GetXsUncertainty(const fastNLO::ENumUncertaintyStyle eNumUnc, bool lNorm = false);
   // std::vector< std::vector<double> > GetXsUncertaintyVec(const fastNLO::ENumUncertaintyStyle eNumUnc, bool lNorm = false, int iprint = 0);
   // void PrintXsUncertaintyVec(fastNLO::ENumUncertaintyStyle, std::string UncName, bool lNorm = false);
   std::vector< std::vector<double> > GetNumUncertaintyVec(const fastNLO::ENumUncertaintyStyle eNumUnc, bool lNorm = false, int iprint = 0);
   void PrintNumUncertaintyVec(fastNLO::ENumUncertaintyStyle, std::string UncName, bool lNorm = false);
   //
   XsUncertainty GetXsUncertainty(const fastNLO::EScaleUncertaintyStyle eScaleUnc, bool lNorm = false, double sclfac = 1.);
   // std::vector< std::vector<double> > GetXsUncertaintyVec(const fastNLO::EScaleUncertaintyStyle eScaleUnc, bool lNorm = false, int iprint = 0, double sclfac = 1.);
   // void PrintXsUncertaintyVec(fastNLO::EScaleUncertaintyStyle, std::string UncName, bool lNorm = false, double sclfac =1.);
   std::vector< std::vector<double> > GetScaleUncertaintyVec(const fastNLO::EScaleUncertaintyStyle eScaleUnc, bool lNorm = false, int iprint = 0, double sclfac = 1.);
   void PrintScaleUncertaintyVec(fastNLO::EScaleUncertaintyStyle, std::string UncName, bool lNorm = false, double sclfac =1.);
   // Specific implementations in fastNLOLHAPDF
   XsUncertainty GetXsUncertainty(const fastNLO::EAsUncertaintyStyle eAsUnc, bool lNorm = false);
   // std::vector< std::vector<double> > GetXsUncertaintyVec(const fastNLO::EAsUncertaintyStyle eAsUnc, bool lNorm = false, int iprint = 0);
   // void PrintXsUncertaintyVec(fastNLO::EAsUncertaintyStyle, std::string UncName, bool lNorm = false);
   std::vector< std::vector<double> > GetAsUncertaintyVec(const fastNLO::EAsUncertaintyStyle eAsUnc, bool lNorm = false, int iprint = 0);
   void PrintAsUncertaintyVec(fastNLO::EAsUncertaintyStyle, std::string UncName, bool lNorm = false);
   //
   XsUncertainty GetXsUncertainty(const fastNLO::EPDFUncertaintyStyle ePDFUnc, bool lNorm = false);
   // std::vector<std::vector<double> > GetXsUncertaintyVec(const fastNLO::EPDFUncertaintyStyle, bool lNorm = false, int iprint = 0);
   // void PrintXsUncertaintyVec(fastNLO::EPDFUncertaintyStyle, std::string UncName, bool lNorm = false);
   std::vector<std::vector<double> > GetPDFUncertaintyVec(const fastNLO::EPDFUncertaintyStyle, bool lNorm = false, int iprint = 0);
   void PrintPDFUncertaintyVec(fastNLO::EPDFUncertaintyStyle, std::string UncName, bool lNorm = false);

   std::vector<LHAPDF::PDFUncertainty>  GetPDFUncertaintyLHAPDF(double cl=100*erf(1/sqrt(2)), bool alternative=false); //!< return PDF uncertainty, formulae taken from LHAPDF6
   std::vector<double> CalcPDFUncertaintyMinus(const std::vector<LHAPDF::PDFUncertainty>& ) const; //!<get vector<double> for PDF-minus uncertainty. Uncertainties are POSITIVE!
   std::vector<double> CalcPDFUncertaintyPlus(const std::vector<LHAPDF::PDFUncertainty>& ) const; //!<get vector<double> for PDF-up uncertainty
   std::vector<double> CalcPDFUncertaintyRelMinus(const std::vector<LHAPDF::PDFUncertainty>& ) const; //!<get vector<double> for relative PDF-minus uncertainty. Uncertainties are NEGATIVE!
   std::vector<double> CalcPDFUncertaintyRelPlus(const std::vector<LHAPDF::PDFUncertainty>& ) const; //!<get vector<double> for relative PDF-up uncertainty
   std::vector<double> CalcPDFUncertaintySymm(const std::vector<LHAPDF::PDFUncertainty>& ) const; //!< get vector<double> for symmetrized PDF uncertainty
   std::vector<double> CalcPDFUncertaintyCentral(const std::vector<LHAPDF::PDFUncertainty>& ) const; //!< get vector<double> for 'new' central value

   // inherited functions
   virtual double EvolveAlphas(double Q) const ;
   virtual bool InitPDF();
   virtual std::vector<double> GetXFX(double xp, double muf) const ;

protected:

   // ---- LHAPDF vars ---- //
   std::string fLHAPDFFilename;
   LHAPDF::PDFSet* PDFSet;
   LHAPDF::PDF* PDF;
   int fnPDFs;
   int fiPDFMember;

   double fchksum;


};

#endif
