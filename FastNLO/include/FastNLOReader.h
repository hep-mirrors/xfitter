
// Author: Daniel Britzger
// DESY, 23/07/2011

//  Version 0.5, 
//
//  History:
//    Version 0, initial version


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


#include <string>
#include <iostream>
#include <cstdio>
#include <vector>
#include "FastNLOBlockB.h"

using namespace std;

class FastNLOReader {

public:
   enum EMuX {
      kMuR			= 0,	// renormalization scale
      kMuF			= 1	// factorization scale
   };
  
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
      kExtern			= 9	// define an external function for your scale
   };

   enum EPDFInterface {
      kLHAPDF			= 0,	// use LHAPDF
      kH1FITTER			= 1	// use H1Fitter for determining the pdflc
   };

   enum EAlphasEvolution {
      kGRV			= 0,
      kNLOJET			= 1,      
      kCTEQpdf			= 2,      
      kLHAPDFInternal		= 4,	// use: double 	LHAPDF::alphasPDF (double Q)
      kQCDNUMInternal		= 5,	// You cannot change alpha_s(Mz) here, but is done within QCDNUM
      kFixed                    = 6     // Always gives back alpha_s(Mz) for testing.
   };

   enum EScaleVariationDefinition {
      kMuRVar			= 0,	// vary only MuR var by a factor 2 up and down
      kMuRMuFSimultaneously	= 1,	// vary MuR and MuF simulataneously up and down
      kMuRTimesMuFVar		= 2,	// vary MuR x MuF up and down by a factor of 2
      kUpDownMax		= 3	// performa a scan for maximum up an down value between 0.5 < cR x cF < 2
   };

   enum EUnits {
      kAbsoluteUnits		= 0,	// calculate the cross section in barn for each publicated bin
      kPublicationUnits		= 1	// calculate the cross section in units as given in the according publication
   };

   enum ESMCalculation {
      kFixedOrder		= 0,	// Fixed Order Calculation
      kThresholdCorrection	= 1,	// Threshold corrections
      kElectroWeakCorrection	= 2,	// Electro weak corrections
      kNonPerturbativeCorrection	= 10	// Hadronisation correction/non-perturbative correction
   };

   static const double TWOPI = 6.28318530717958647692528;
   static const double TWOPISQR = 39.47841760435743447533796;

protected:

   static const int tablemagicno	= 1234567890;
   string ffilename;

   // ---- scale variations ---- //
   int fScalevar;
   double fScaleFacMuR;
   double fScaleFacMuF;

   // ---- LHAPDF vars ---- //
   string fLHAPDFfilename;
   //string fLHAPDFpath;
   int fnPDFs;
   int fiPDFSet;

   EPDFInterface	fPDFInterface;
   EAlphasEvolution	fAlphasEvolution;
   EScaleVariationDefinition fScaleVariationDefinition;
   EScaleFunctionalForm fMuRFunc;
   EScaleFunctionalForm fMuFFunc;
   EUnits		fUnits;
   double (*Fct_MuR)(double,double);			// Function, if you define your functional form for your scale external
   double (*Fct_MuF)(double,double);			// Function, if you define your functional form for your scale external
   //    ECalculationOrder	fOrder;
   // ---- alpha_s vars ---- //
   double fAlphasMz;

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
   vector < vector < FastNLOBlockB* > > BBlocksSMCalc;	// BlockB's for SM corrections (IContrFlag1 = 2) [Model(i~ContrFlag2)][contribution]
							//  e.g. model = th. corr LO and NLO, e/w LO and NLO and NNLO, 0->fixed order
							//  what to do with each contribution is defined in IContrFlag3
   vector < vector < FastNLOBlockB* > > BBlocksNewPhys;		// BlockB's for New physics corrections [model][contribution]

   vector < vector < bool > > bUseSMCalc;		// switch calclations ON/OFF
   vector < vector < bool > > bUseNewPhys;		// switch calclations ON/OFF

   // ---- Cross sections ---- //
   vector < double > XSection_LO;
   vector < double > XSection;
   // k-factor
   vector < double > kFactor;

   // ----  reference tables ---- //
   // v2.0
   vector < double > XSectionRef;
   // v2.0+ MuVar
   //vector < double > XSectionMuVarRef;
   vector < double > XSectionRefMixed;
   vector < double > XSectionRef_s1;
   vector < double > XSectionRef_s2;

private:

   void Init() ;
   void ReadTable();
   void StripWhitespace(string* s);

   void ReadBlockA1(istream *table);
   void ReadBlockA2(istream *table);
   void ReadBlockB(istream *table);
  
   void PrintBlockA1();
   void PrintBlockA2();

   void InitLHAPDF();
   void FillBlockBPDFLCsDISv20( FastNLOBlockB* B );
   void FillBlockBPDFLCsDISv21( FastNLOBlockB* B );
   void FillBlockBPDFLCsHHCv20( FastNLOBlockB* B );
   void FillBlockBPDFLCsHHCv21( FastNLOBlockB* B );
   void CalcAposterioriScaleVariation();
   vector<double> GetXFX(double x, double muf);
   vector<double> CalcPDFLinearCombDIS(vector<double> pdfx1, int NSubproc );
   vector<double> CalcPDFLinearCombHHC(vector<double> pdfx1, vector<double> pdfx2, int NSubproc );
   void FillAlphasCacheInBlockBv20( FastNLOBlockB* B );
   void FillAlphasCacheInBlockBv21( FastNLOBlockB* B );
   double GetAlphas(double Q);

   void CalcReferenceCrossSection();

   double GetAlphasLHAPDF(double Q);
   double GetAlphasQCDNUM(double Q);
   double GetAlphasNLOJET(double Q, double alphasMz);
   double GetAlphasGRV(double Q, double alphasMz);
   double GetAlphasNewGRV(double Q, double alphasMz);
   double GetAlphasCTEQpdf(double Q, double alphasMz);
   double GetAlphasFixed(double Q, double alphasMz);
  
   double CalcMu(FastNLOReader::EMuX kMuX, double scale1 , double scale2 , double scalefactor);
   double FuncMixedOver1 ( double scale1 , double scale2 ) ;
   double FuncMixedOver2 ( double scale1 , double scale2 ) ;
   double FuncMixedOver4 ( double scale1 , double scale2 ) ;
   double FuncLinearMean ( double scale1 , double scale2 ) ;
   double FuncLinearSum ( double scale1 , double scale2 ) ;
   double FuncMax ( double scale1 , double scale2 ) ;
   double FuncMin ( double scale1 , double scale2 ) ;

   void CalcCrossSectionv21(FastNLOBlockB* B , bool IsLO = false );
   void CalcCrossSectionv20(FastNLOBlockB* B , bool IsLO = false);
 
public:
   FastNLOReader(void);
   FastNLOReader(string filename);
   ~FastNLOReader(void);

   void SetFilename(string filename) ;
   void InitScalevariation();
   void SetLHAPDFfilename( string filename ) { fLHAPDFfilename = filename; };
   //void SetLHAPDFpath( string path ) { fLHAPDFpath = path; };
   void SetLHAPDFset( int set ) { fiPDFSet = set; };
   void SetAlphasMz( double AlphasMz , bool ReCalcCrossSection = false );
   void SetPDFInterface( EPDFInterface PDFInterface)	{ fPDFInterface = PDFInterface; };
   void SetAlphasEvolution( EAlphasEvolution AlphasEvolution );
   void SetScaleVariationDefinition( EScaleVariationDefinition ScaleVariationDefinition ) { fScaleVariationDefinition = ScaleVariationDefinition ; cout << "not implemented yet."<<endl;}; // no impact yet.
   void SetUnits( EUnits Unit );
   //void SetCalculationOrder( ECalculationOrder order ){ fOrder = order;};
   void SetContributionON( ESMCalculation eCalc , unsigned int Id , bool SetOn = true );	// Set contribution On/Off. Look for Id of this contribution during initialization.

   // ---- setters for scales of MuVar tables ---- //
   void SetMuRFunctionalForm( EScaleFunctionalForm func , bool ReFillCache = true );	// Set the functional form of Mu_R
   void SetMuFFunctionalForm( EScaleFunctionalForm func , bool ReFillCache = true );	// Set the functional form of Mu_F
   void SetFunctionalForm( EScaleFunctionalForm func , FastNLOReader::EMuX kMuX );	// Set functional form of MuX
   void SetScaleFactorMuR( double fac , bool ReFillCache = true );			// Set scale factor for MuR
   void SetScaleFactorMuF( double fac , bool ReFillCache = true );			// Set scale factor for MuF
   void SetExternalFuncForMuR( double (*Func)(double,double) , bool ReFillCache = true );	// Set external function for scale calculation (optional)
   void SetExternalFuncForMuF( double (*Func)(double,double) , bool ReFillCache = true );	// Set external function for scale calculation (optional)


   // ---- setters for scale variation in v2.0 tables  ---- //
   double SetScaleVariation(int scalevar , bool ReFillCache = true);			// choose the scale variation table
  
   // ---- Pdf interface ---- //
   void FillPDFCache( bool ReCalcCrossSection = false );				// Prepare for recalculation of cross section with 'new'/updated pdf.

   // ---- alphas cache ---- //
   void FillAlphasCache();								// prepare for recalculation of cross section with new alpha_s value.

   // ---- Getters ---- //
   vector < double > GetCrossSection();
   vector < double > GetReferenceCrossSection();
   vector < double > GetKFactors();

   int GetNcontrib() { return Ncontrib; };
   int GetIExpUnit() { return Ipublunits; };				// exponent of xs units (like -12 for pb)
   string GetScenarioName() { return ScenName; };			// Get Scenario/Table name
   vector < string > GetScenarioDescription() { return ScDescript; };	// Get Description of scenario
   double GetCMSEnergy() { return Ecms; };				// Get center of mass energy
   int GetILOord() { return ILOord; };					// Get number of alpha_s in leading order (1 or 2 usually)
   int GetNObsBins() { return NObsBin; };				// Get number of measured bins
   int GetNDiffBin() { return NDim; };					// Get number of differential measurement. 1: single differential; 2: double differential
   vector < int > GetRapidityIndex() { return RapIndex;};		// Get rapidity indices
   vector < string > GetDimensionLabel() { return DimLabel;};		// Get label for measured dimensions
   vector < int > GetIDiffBin() { return IDiffBin;};			// Get number of differential bins
   vector < vector < double > > GetLowBinEdge() { return LoBin; };	// Get Lower Bin edge [ObsBin][DiffBin]
   vector < vector < double > > GetUpBinEdge() { return UpBin; };	// Get Upper Bin edge [ObsBin][DiffBin]
   vector < double > GetBinSize() { return BinSize; };			// Get Binsize = BinSizeDim1 < * BinSizeDim2 >
   int IsNormalized() { return INormFlag; };				// Normalized Cross sections?
   //   string DenomTable;
   //   vector <int> IDivLoPointer;
   //   vector <int> IDivUpPointer;

   // Getters about scale-interpolations
   string GetScaleDescription() { return BBlocksSMCalc[0][1]->ScaleDescript[0][0]; };		// Description of renormalization and facorization scale choice
   int GetNScaleVariations();									// Get number of available scale variations
   vector < double > GetScaleFactors();								// Get list of available scale factors
  

   void CalcCrossSection();
   void PrintTableInfo(const int iprint = 0);
   void PrintDataCrossSections();
   void PrintFastNLOTableConstants(const int iprint = 2);
   void PrintCrossSections();
   void PrintCrossSectionsLikeFreader();
   void PrintCrossSectionsWithReference();

   static const string fOrdName[4];
   static const string fCorrName[11];
   static const string fNPName[10];
   static const string fNSDep[4];

private:
   static int WelcomeOnce;

};

#endif
