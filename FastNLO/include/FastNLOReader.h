// Author: Daniel Britzger
// DESY, 23/07/2011

//  Version 0.1, 
//
//  History:
//    Version 0, initial version


#ifndef FASTNLOREADER
#define FASTNLOREADER

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                                                                      //
//  FastNLOReader is a standalone code for reading                      //
//  FastNLO tables of version 2.0 for DIS, pp, and ppbar                //
//  processes. It is also optimized for an integration into             //
//  the H1Fitter project.                                               //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <string>
#include <iostream>
#include <cstdio>
#include <vector>
#include <FastNLOBlockB.h>

using namespace std;

class FastNLOReader { //: public TObject {
 private:
   void InitFastNLOREader(void);     // initialize all data members
 public:
   enum EPDFInterface {
      kLHAPDF	= 0,	// use LHAPDF
      kH1FITTER = 1	// use H1Fitter for determining the pdflc
   };

   enum EAlphasEvolution {
      kGRV = 0,
      kNLOJET = 1,      
      kCTEQpdf = 2,      
      kFastNLO = 3,
      kLHAPDFInternal = 4		// use: double 	LHAPDF::alphasPDF (double Q)
   };

   enum EScaleVariationDefinition {
      kMuRVar		= 0,	// vary only MuR var by a factor 2 up and down
      kMuRMuFSimultaneously	= 1,	// vary MuR and MuF simulataneously up and down
      kMuRTimesMuFVar	= 2,	// vary MuR x MuF up and down by a factor of 2
      kUpDownMax	= 3	// performa a scan for maximum up an down value between 0.5 < cR x cF < 2
   };

  static const double TWOPI = 6.28318530717958647692528;
  static const double TWOPISQR = 39.47841760435743447533796;

 protected:
  // put all the stuff in here!

  static const int tablemagicno	= 1234567890;
  string ffilename;
  int fScalevar;

  // ---- LHAPDF vars ---- //
  string fLHAPDFfilename;
  string fLHAPDFpath;
  int fnPDFs;
  int fiPDFSet;

  EPDFInterface	fPDFInterface;
  EAlphasEvolution fAlphasEvolution;
  EScaleVariationDefinition fScaleVariationDefinition;

  // ---- LHAPDF vars ---- //
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
  int NScDescript;
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
  FastNLOBlockB* BlockB_LO;
  FastNLOBlockB* BlockB_NLO;
  FastNLOBlockB* BlockB_LO_Ref;
  FastNLOBlockB* BlockB_NLO_Ref;

  // ---- Cross sections ---- //
  // v2.0
  vector < double > XSection_LO;
  vector < double > XSection;
  // v2.0 with two scales
  vector < double > XSection2Scales_LO;
  vector < double > XSection2Scales;
  // v2.0+ MuVar
  vector < double > XSectionMuVar_LO;
  vector < double > XSectionMuVar;

  // ----  reference tables ---- //
  // v2.0
  vector < double > XSectionRef;
  // v2.0 with two scales
  vector < double > XSection2ScalesRef;
  // v2.0+ MuVar
  //vector < double > XSectionMuVarRef;
  vector < double > XSectionRefMixed;
  vector < double > XSectionRefQ2;
  vector < double > XSectionRefMufQ2MuRMixed;

  // ----  reference tables ---- //

 private:

  void ReadTable();

  void ReadBlockA1(istream *table);
  void ReadBlockA2(istream *table);
  void ReadBlockB(istream *table);
  
  void PrintBlockA1();
  void PrintBlockA2();

  void InitLHAPDF();
  void FillBlockBPDFLCsWithLHAPDF( FastNLOBlockB* B );
  void FillBlockBPDFLCsWithH1Fitter( FastNLOBlockB* B );
  vector<double> CalcPDFLinearComb(vector<double> pdfx1, vector<double> pdfx2, int IPDFdef1, int IPDFdef2, int NSubproc );
  vector<double> CalcPDFLinearCombPPMuVar(vector<double> pdfx1, vector<double> pdfx2 );
  vector<double> CalcPDFLinearCombPPbarMuVar(vector<double> pdfx1, vector<double> pdfx2 );
  vector<double> CalcPDFLinearCombDIS(vector<double> pdfx1, int NSubproc );
  void FillAlphasCacheInBlockB( FastNLOBlockB* B );
  double GetAlphas(double Q);

  void CalcReferenceCrossSection();

  double GetAlphasLHAPDF(double Q);
  double GetAlphasNLOJET(double Q, double alphasMz);
  double GetAlphasGRV(double Q, double alphasMz);
  double GetAlphasCTEQpdf(double Q, double alphasMz);
  double GetAlphasFastNLO(double Q, double alphasMz);
  
  
 protected:
//    TUnfold(void);              // for derived classes
//    virtual Double_t DoUnfold(void);     // the unfolding algorithm
//    virtual void ClearResults(void);     // clear all results
  
public:
  FastNLOReader(void);
  FastNLOReader(string filename);

  void SetFilename(string filename) { ffilename = filename;};
  void SetLHAPDFfilename( string filename ) { fLHAPDFfilename = filename; };
  void SetLHAPDFpath( string path ) { fLHAPDFpath = path; };
  void SetLHAPDFset( int set ) { fiPDFSet = set; };
  void SetAlphasMz( double AlphasMz );
  void SetPDFInterface( EPDFInterface PDFInterface)	{ fPDFInterface = PDFInterface; };
  void SetAlphasEvolution( EAlphasEvolution AlphasEvolution ) { fAlphasEvolution = AlphasEvolution; if (AlphasEvolution==kLHAPDFInternal) cout << "Warning. You cannot change the Alpha_s value."<<endl; };
  void SetScaleVariationDefinition( EScaleVariationDefinition ScaleVariationDefinition ) { fScaleVariationDefinition = ScaleVariationDefinition ;};

  double SetScaleVariation(int scalevar);
  void FillPDFCache();
  void FillAlphasCache();

  vector < double > GetXSection();
  vector < double > GetReferenceXSection();

  void CalcCrossSection();
  void Print();

  //  ClassDef(FastNLOReader, 0) // FastNLOReader. Standalone code for reading FastNLO v2.0 and v2.0+ tables

};

#endif
