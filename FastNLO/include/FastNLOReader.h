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
    kQuadraticSumOver4		= 4,	// e.g. mu^2 = ( Q^2 + pt^2 ) / 4 
    kLinearMean			= 5,	// e.g. mu^2 = (( Q + pt ) / 2 )^2
    kLinearSum			= 6,	// e.g. mu^2 = (( Q + pt ))^2
    kScaleMax			= 7,	// e.g. mu^2 = max( Q^2, pt^2)
    kScaleMin			= 8,	// e.g. mu^2 = min( Q^2, pt^2) 
    kExpProd2                   = 9,    // e.g. mu^2 = (scale1 * exp(0.3 * scale2)) ^2
    kExtern			= 10	// define an external function for your scale
  };

  enum EPDFInterface {
    kLHAPDF			= 0,	// use LHAPDF
    kQCDNUM			= 1,	// use QCDNUM for the pdf grid
    kH1Fitter			= 2,	// use H1Fitter/HeraFitter for the pdf grid
    kDiffPDF			= 3	// use DiffPDF
  };

  enum EAlphasEvolution {
    kGRV			= 0,
    kNLOJET			= 1,      
    kCTEQpdf			= 2,      
    kLHAPDFAs			= 4,	// use: double 	LHAPDF::alphasPDF (double Q)
    kQCDNUMAs			= 5,	// You cannot change alpha_s(Mz) here, but is done within QCDNUM
    kH1FitterAs			= 6,	// Alpha_s routine from H1Fitter
    kFixed			= 7     // Always gives back alpha_s(Mz) for testing.
  };

  enum EScaleVariationDefinition {
    kMuRVar			= 0,	// vary only MuR var by a factor 2 up and down
    kMuRMuFSimultaneously	= 1,	// vary MuR and MuF simulataneously up and down
    kMuRTimesMuFVar		= 2,	// vary MuR x MuF up and down by a factor of 2
    kUpDownMax			= 3	// perform a scan for maximum up an down value between 0.5 < cR x cF < 2
  };

  enum EUnits {
    kAbsoluteUnits		= 0,	// calculate the cross section in barn for each publicated bin
    kPublicationUnits		= 1	// calculate the cross section in units as given in the according publication
  };

  // Corresponds to IContrFlag1 in v2 table definition
  enum ESMCalculation {
    kFixedOrder		        = 0,	// Fixed order calculation (pQCD)
    kThresholdCorrection	= 1,	// Threshold corrections
    kElectroWeakCorrection	= 2,	// Electroweak corrections
    kNonPerturbativeCorrection	= 3	// Non-perturbative corrections|Hadronisation corrections
  };
  
  // Corresponds to IContrFlag2 in v2 table definition
  enum ESMOrder {
    kLeading		        = 0,	// LO,   1-loop, LO MC
    kNextToLeading        	= 1,	// NLO,  2-loop, NLO MC
    kNextToNextToLeading	= 2	// NNLO, 3-loop, NNLO MC
  };
  
  static const double TWOPI = 6.28318530717958647692528;
  static const double TWOPISQR = 39.47841760435743447533796;

  typedef double(*mu_func)(double,double);

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

   // ---- vars for diffractive tables ---- //
   double fxpom;
   double fzmin;
   double fzmax;

  EPDFInterface	fPDFInterface;
  EAlphasEvolution	fAlphasEvolution;
  EScaleVariationDefinition fScaleVariationDefinition;
  EScaleFunctionalForm fMuRFunc;
  EScaleFunctionalForm fMuFFunc;
  EUnits		fUnits;
  mu_func Fct_MuR;				// Function, if you define your functional form for your scale external
  mu_func Fct_MuF;				// Function, if you define your functional form for your scale external

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
  vector < vector < FastNLOBlockB* > > BBlocksSMCalc;	// BlockB's for SM corrections
  vector < vector < FastNLOBlockB* > > BBlocksNewPhys;	// BlockB's for New physics corrections

  vector < vector < bool > > bUseSMCalc;		// switch calculations ON/OFF
  vector < vector < bool > > bUseNewPhys;		// switch calculations ON/OFF

  // ---- Cross sections ---- //
  vector < double > XSection_LO;
  vector < double > XSection;
  vector < double > kFactor;

  // ----  reference tables ---- //
  vector < double > XSectionRef;
  vector < double > XSectionRefMixed;
  vector < double > XSectionRef_s1;
  vector < double > XSectionRef_s2;

 
public:

  FastNLOReader(void);
  FastNLOReader(string filename);
  ~FastNLOReader(void);

  void SetFilename(string filename) ;
  void InitScalevariation();
  void SetLHAPDFfilename( string filename );
  //void SetLHAPDFpath( string path ) { fLHAPDFpath = path; };
  void SetLHAPDFset( int set );
  void PrintCurrentLHAPDFInformation() const;
  void SetAlphasMz( double AlphasMz , bool ReCalcCrossSection = false );
  void SetPDFInterface( EPDFInterface PDFInterface)	{ fPDFInterface = PDFInterface; };
  void SetAlphasEvolution( EAlphasEvolution AlphasEvolution );
  void SetScaleVariationDefinition( EScaleVariationDefinition ScaleVariationDefinition ) { fScaleVariationDefinition = ScaleVariationDefinition ; cout << "not implemented yet."<<endl;}; // no impact yet.
  void SetUnits( EUnits Unit );
  //void SetCalculationOrder( ECalculationOrder order ){ fOrder = order;};
  void SetContributionON( ESMCalculation eCalc , unsigned int Id , bool SetOn = true, bool Verbose = false );	// Set contribution On/Off. Look for Id of this contribution during initialization.
  int ContrId( const ESMCalculation eCalc, const ESMOrder eOrder ) const;
  void SetGRVtoPDG2011_2loop(bool print);

  // ---- setters for scales of MuVar tables ---- //
  void SetMuRFunctionalForm( EScaleFunctionalForm func , bool ReFillCache = true );	// Set the functional form of Mu_R
  void SetMuFFunctionalForm( EScaleFunctionalForm func , bool ReFillCache = true );	// Set the functional form of Mu_F
  void SetFunctionalForm( EScaleFunctionalForm func , FastNLOReader::EMuX kMuX );	// Set functional form of MuX
  void SetScaleFactorMuR( double fac , bool ReFillCache = true );			// Set scale factor for MuR
  void SetScaleFactorMuF( double fac , bool ReFillCache = true );			// Set scale factor for MuF
  void SetExternalFuncForMuR( mu_func , bool ReFillCache = true );			// Set external function for scale calculation (optional)
  void SetExternalFuncForMuF( mu_func , bool ReFillCache = true );			// Set external function for scale calculation (optional)


  // ---- setters for scale variation in v2.0 tables  ---- //
  double SetScaleVariation(int scalevar , bool ReFillCache = true);			// choose the scale variation table
  

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
  EScaleFunctionalForm GetMuRFunctionalForm() const { return fMuRFunc; };
  EScaleFunctionalForm GetMuFFunctionalForm() const { return fMuFFunc; };
  EPDFInterface GetPDFInterface() const  { return fPDFInterface; };
  EAlphasEvolution GetAlphasEvolution() const { return fAlphasEvolution; };
  EScaleVariationDefinition GetScaleVariationDefinition() const { return fScaleVariationDefinition; };
  EUnits GetUnits() const{ return fUnits; };
  mu_func GetExternalFuncForMuR(){ return Fct_MuR; };
  mu_func GetExternalFuncForMuF(){ return Fct_MuF; };
  double GetAlphasMz() const { return fAlphasMz; };
  double GetScaleFactorMuR() const { return fScaleFacMuR;};
  double GetScaleFactorMuF() const { return fScaleFacMuF;};
  int GetScaleVariation() const { return fScalevar; };
  int GetIPDFSet() const {return fiPDFSet;};
  int GetNPDFSets() const {return fnPDFs;};


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
  string GetScaleDescription() const { return BBlocksSMCalc[0][1]->ScaleDescript[0][0]; };		// Description of renormalization and facorization scale choice
  int GetNScaleVariations() const;					// Get number of available scale variations
  vector < double > GetScaleFactors() const;				// Get list of available scale factors
  

  // ---- Print outs ---- //
  void PrintTableInfo(const int iprint = 0) const;
  void PrintFastNLOTableConstants(const int iprint = 2) const;
  void PrintCrossSections() const ;
  void PrintCrossSectionsDefault(vector<double> kthc = vector<double>() ) const ;
  void PrintCrossSectionsWithReference();
  void PrintCrossSectionsData() const;
  void PrintFastNLODemo();


  // ---- human readable strings ---- //
  static const string fContrName[20];
  static const string fOrdName[4][4];
  static const string fNSDep[4];


protected:

  void Init() ;
  void ReadTable();
  void InitMembers();
  void StripWhitespace(string* s);

  void ReadBlockA1(istream *table);
  void ReadBlockA2(istream *table);
  void ReadBlockB(istream *table);

  void PrintBlockA1() const;
  void PrintBlockA2() const;

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
  double CalcAlphas(double Q);

  void CalcReferenceCrossSection();

  double CalcAlphasLHAPDF(double Q);
  double CalcAlphasQCDNUM(double Q);
  double CalcAlphasNLOJET(double Q, double alphasMz);
  double CalcAlphasGRV(double Q, double alphasMz);
  double CalcAlphasNewGRV(double Q, double alphasMz);
  double CalcAlphasCTEQpdf(double Q, double alphasMz);
  double CalcAlphasFixed(double Q, double alphasMz);
  
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

protected:
  static int WelcomeOnce;

};

#endif
