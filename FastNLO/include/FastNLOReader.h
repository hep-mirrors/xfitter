// Author: Daniel Britzger
// DESY, 23/07/2011

#ifndef FASTNLOREADER
#define FASTNLOREADER


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  FastNLOReader                                                       //
//                                                                      //
//  FastNLOReader is a standalone code for reading                      //
//  FastNLO tables of version 2.0 for DIS processes                     //
//  It is also optimized for an integration into                        //
//  the H1Fitter project.                                               //
//                                                                      //
//  FastNLO is developed by                                             //
//    D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch         //
//    (publication in preparation)                                      //
//    http://projects.hepforge.org/fastnlo                              //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "speaker.h"
#include <string>
#include <vector>
#include "FastNLOBlockB.h"

namespace fastNLO {
   enum EScaleFunctionalForm {
      kScale1			= 0,	// e.g. mu^2 = Q^2 
      kScale2			= 1,	// e.g. mu^2 = pt^2 
      kQuadraticSum		= 2,	// e.g. mu^2 = ( Q^2 + pt^2 )
      kQuadraticMean		= 3,	// e.g. mu^2 = ( Q^2 + pt^2 ) / 2 
      kQuadraticSumOver4	= 4,	// e.g. mu^2 = ( Q^2 + pt^2 ) / 4 
      kLinearMean		= 5,	// e.g. mu^2 = (( Q + pt ) / 2 )^2
      kLinearSum		= 6,	// e.g. mu^2 = (( Q + pt ))^2
      kScaleMax			= 7,	// e.g. mu^2 = max( Q^2, pt^2)
      kScaleMin			= 8,	// e.g. mu^2 = min( Q^2, pt^2) 
      kExpProd2			= 9,    // e.g. mu^2 = (scale1 * exp(0.3 * scale2)) ^2
      kExtern			= 10	// define an external function for your scale
   };
   enum ESMCalculation {
      kFixedOrder		= 0,	// Fixed order calculation (pQCD)
      kThresholdCorrection	= 1,	// Threshold corrections
      kElectroWeakCorrection	= 2,	// Electroweak corrections
      kNonPerturbativeCorrection= 3	// Non-perturbative corrections|Hadronisation corrections
   };
   enum ESMOrder {
      kLeading		        = 0,	// LO,   1-loop, LO MC
      kNextToLeading        	= 1,	// NLO,  2-loop, NLO MC
      kNextToNextToLeading	= 2	// NNLO, 3-loop, NNLO MC
   };
   enum EUnits {
      kAbsoluteUnits		= 0,	// calculate the cross section in barn for each publicated bin
      kPublicationUnits		= 1	// calculate the cross section in units as given in the according publication
   };
}

using namespace std;
using namespace fastNLO;

class FastNLOReader : public PrimalScream {

public:

  typedef double(*mu_func)(double,double);
  enum EMuX {
    kMuR			= 0,	// renormalization scale
    kMuF			= 1	// factorization scale
  };
  static const double TWOPI = 6.28318530717958647692528;
  static const double TWOPISQR = 39.47841760435743447533796;

public:

  FastNLOReader(string filename);
  FastNLOReader(const FastNLOReader& fnlo);
  virtual ~FastNLOReader(void);

  void SetFilename(string filename) ;
  void InitScalevariation();
  void SetUnits( fastNLO::EUnits Unit );
  void SetContributionON( fastNLO::ESMCalculation eCalc , unsigned int Id , bool SetOn = true);	// Set contribution On/Off. Look for Id of this contribution during initialization.
  int ContrId( const fastNLO::ESMCalculation eCalc, const fastNLO::ESMOrder eOrder ) const;

  // ---- setters for scales of MuVar tables ---- //
  void SetMuRFunctionalForm( fastNLO::EScaleFunctionalForm func);// Set the functional form of Mu_R
  void SetMuFFunctionalForm( fastNLO::EScaleFunctionalForm func , bool ReFillCache = true);// Set the functional form of Mu_F
  void SetFunctionalForm( fastNLO::EScaleFunctionalForm func , FastNLOReader::EMuX kMuX);// Set functional form of MuX
  bool SetScaleFactorsMuRMuF( double xmur, double xmuf, bool ReFillCache = true);// Set scale factors for MuR and MuF
  void SetExternalFuncForMuR( mu_func);						// Set external function for scale calculation (optional)
  void SetExternalFuncForMuF( mu_func , bool ReFillCache = true);		// Set external function for scale calculation (optional)


  // ---- Pdf interface ---- //
  void FillPDFCache( bool ReCalcCrossSection = false );					// Prepare for recalculation of cross section with 'new'/updated pdf.

   // ---- alphas cache ---- //
   void FillAlphasCache();								// prepare for recalculation of cross section with new alpha_s value.

  // ---- Do the cross section calculation ---- //
  void CalcCrossSection();


  // ---- Getters for results---- //
  vector < double > GetCrossSection();
  vector < double > GetReferenceCrossSection();
  vector < double > GetKFactors();


  // ---- Getters for FastNLOReader member variables ---- //
  fastNLO::EScaleFunctionalForm GetMuRFunctionalForm() const { return fMuRFunc; };
  fastNLO::EScaleFunctionalForm GetMuFFunctionalForm() const { return fMuFFunc; };
  fastNLO::EUnits GetUnits() const{ return fUnits; };
  mu_func GetExternalFuncForMuR(){ return Fct_MuR; };
  mu_func GetExternalFuncForMuF(){ return Fct_MuF; };
  double GetScaleFactorMuR() const { return fScaleFacMuR;};
  double GetScaleFactorMuF() const { return fScaleFacMuF;};
  int GetScaleVariation() const { return fScalevar; };

  // ---- Getters for FastNLO table constants ---- //
  int GetNcontrib() const { return Ncontrib; };
  int GetIExpUnit() const { return Ipublunits; };			// exponent of xs units (like -12 for pb)
  string GetScenarioName() const { return ScenName; };			// Get Scenario/Table name
  vector < string > GetScenarioDescription() const { return ScDescript; };	// Get Description of scenario
  double GetCMSEnergy() const { return Ecms; };				// Get center of mass energy
  int GetILOord() const { return ILOord; };				// Get number of alpha_s in leading order (1 or 2 usually)
  int GetNObsBins() const { return NObsBin; };				// Get number of measured bins
  int GetNDiffBin() const { return NDim; };				// Get number of differential measurement. 1: single differential; 2: double differential
  vector < int > GetRapidityIndex() const { return RapIndex;};		// Get rapidity indices
  vector < string > GetDimensionLabel() const { return DimLabel;};	// Get label for measured dimensions
  vector < int > GetIDiffBin() const { return IDiffBin;};		// Get number of differential bins
  vector < vector < double > > GetLowBinEdge() const { return LoBin; };	// Get Lower Bin edge [ObsBin][DiffBin]
  vector < vector < double > > GetUpBinEdge() const { return UpBin; };	// Get Upper Bin edge [ObsBin][DiffBin]
  vector < double > GetBinSize() const { return BinSize; };		// Get Binsize = BinSizeDim1 < * BinSizeDim2 >
  int IsNormalized() const { return INormFlag; };			// Normalized Cross sections?
  string GetScaleDescription(int scalen=0) const { return BBlocksSMCalc[0][0]->ScaleDescript[0][scalen]; };		// Description of renormalization and facorization scale choice
  int GetNScaleVariations() const;					// Get number of available scale variations
  vector < double > GetScaleFactors() const;				// Get list of available scale factors
   bool GetIsFlexibleScaleTable() const { return BBlocksSMCalc[0][0]->NScaleDep >= 3; } // Get, if this table is a 'flexible scale' table or not.


  // ---- Print outs ---- //
  void PrintTableInfo(const int iprint = 0) const;
  void PrintFastNLOTableConstants(const int iprint = 2) const;
  void PrintCrossSections() const ;
  void PrintCrossSectionsDefault(vector<double> kthc = vector<double>() ) const ;
  void PrintCrossSectionsWithReference();
  void PrintCrossSectionsData() const;
  void PrintFastNLODemo();


protected:

  void Init() ;
  void ReadTable();
  void StripWhitespace(string* s);

  void ReadBlockA1(istream *table);
  void ReadBlockA2(istream *table);
  void ReadBlockB(istream *table);

  void PrintBlockA1() const;
  void PrintBlockA2() const;

  void PrintScaleSettings(EMuX kMuX=kMuR);
  void FillBlockBPDFLCsDISv20( FastNLOBlockB* B );
  void FillBlockBPDFLCsDISv21( FastNLOBlockB* B );
  void FillBlockBPDFLCsHHCv20( FastNLOBlockB* B );
  void FillBlockBPDFLCsHHCv21( FastNLOBlockB* B );
  void CalcAposterioriScaleVariation();
  vector<double> CalcPDFLinearCombDIS(vector<double> pdfx1, int NSubproc );
  vector<double> CalcPDFLinearCombHHC(vector<double> pdfx1, vector<double> pdfx2, int NSubproc );
  void FillAlphasCacheInBlockBv20( FastNLOBlockB* B );
  void FillAlphasCacheInBlockBv21( FastNLOBlockB* B );
  double CalcAlphas(double Q);

  void CalcReferenceCrossSection();
  
  double CalcMu(FastNLOReader::EMuX kMuX, double scale1 , double scale2 , double scalefactor);
  double FuncMixedOver1 ( double scale1 , double scale2 ) ;
  double FuncMixedOver2 ( double scale1 , double scale2 ) ;
  double FuncMixedOver4 ( double scale1 , double scale2 ) ;
  double FuncLinearMean ( double scale1 , double scale2 ) ;
  double FuncLinearSum ( double scale1 , double scale2 ) ;
  double FuncMax ( double scale1 , double scale2 ) ;
  double FuncMin ( double scale1 , double scale2 ) ;
  double FuncExpProd2 ( double scale1 , double scale2 ) ;

  void CalcCrossSectionv21(FastNLOBlockB* B , bool IsLO = false );
  void CalcCrossSectionv20(FastNLOBlockB* B , bool IsLO = false);
   
   FastNLOBlockB* B_NLO() { return BBlocksSMCalc[0][1]; };
   FastNLOBlockB* B_LO() { return BBlocksSMCalc[0][0]; };
   FastNLOBlockB* B_ThC(int n=0) { 
      if ( BBlocksSMCalc[fastNLO::kThresholdCorrection].empty() ) return NULL;
      else return BBlocksSMCalc[fastNLO::kThresholdCorrection][n]; };

   // virtual functions for the user interface
   virtual void InitPDF() = 0;
   virtual vector<double> GetXFX(double x, double muf ) const = 0;
   virtual double EvolveAlphas(double Q) const = 0;

   // ---- setters for scale variation in v2.0 tables  ---- //
   double SetScaleVariation( int scalevar , bool ReFillCache = true , bool FirstCall=false);// Choose the MuF scale variation table

   // ---- human readable strings ---- //
   static const string fContrName[20];
   static const string fOrdName[4][4];
   static const string fNSDep[6];

protected:

  static const int tablemagicno	= 1234567890;
  static int WelcomeOnce;
  string ffilename;
  int fScalevar;
  double fScaleFacMuR;
  double fScaleFacMuF;
  fastNLO::EScaleFunctionalForm fMuRFunc;
  fastNLO::EScaleFunctionalForm fMuFFunc;
  fastNLO::EUnits		fUnits;
  mu_func Fct_MuR;				// Function, if you define your functional form for your scale external
  mu_func Fct_MuF;				// Function, if you define your functional form for your scale external
  vector < vector < bool > > bUseSMCalc;		// switch calculations ON/OFF
  vector < vector < bool > > bUseNewPhys;		// switch calculations ON/OFF

  // ---- Block A1 ---- //
  int Itabversion;
  string ScenName;
  int Ncontrib;
  int Nmult;
  int Ndata;
  int NuserString;
  int NuserInt;
  int NuserFloat;
  int Imachine;

  // ---- Block A2 ---- //
  int Ipublunits;
  vector < int > bla;
  vector <string> ScDescript;
  double Ecms;
  int ILOord;
  int NObsBin;
  int NDim;
  vector <int> RapIndex;
  vector <string> DimLabel;
  vector <int> IDiffBin;
  vector < vector <double> > LoBin;
  vector < vector <double> > UpBin;
  vector <double> BinSize;
  int INormFlag;
  string DenomTable;
  vector <int> IDivLoPointer;
  vector <int> IDivUpPointer;

  // ---- Block B ---- //
  FastNLOBlockB* BlockB_Data;
  FastNLOBlockB* BlockB_LO_Ref;
  FastNLOBlockB* BlockB_NLO_Ref;
  vector < vector < FastNLOBlockB* > > BBlocksSMCalc;	// BlockB's for SM corrections
  vector < vector < FastNLOBlockB* > > BBlocksNewPhys;	// BlockB's for New physics corrections

  // ---- Cross sections ---- //
  vector < double > XSection_LO;
  vector < double > XSection;
  vector < double > kFactor;

  // ----  reference tables ---- //
  vector < double > XSectionRef;
  vector < double > XSectionRefMixed;
  vector < double > XSectionRef_s1;
  vector < double > XSectionRef_s2;

};


#endif
