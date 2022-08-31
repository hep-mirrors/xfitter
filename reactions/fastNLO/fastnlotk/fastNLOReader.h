#ifndef __fastNLOReader__
#define __fastNLOReader__

#include "fastNLOTable.h"
#include "fastNLOPDFLinearCombinations.h"

// ---- Getters for results---- //
struct XsUncertainty {
   //! Struct for returning vectors with cross section and relative uncertainty
   // keep definition of this class outside of fastNLOReader, because of python wrapper
   std::vector < double > xs;
   std::vector < double > dxsl;
   std::vector < double > dxsu;
};


class fastNLOReader : public fastNLOTable , public fastNLOPDFLinearCombinations {
   //!
   //! fastNLOReader.
   //! Abstract base class for evaluation of fastNLO tables.
   //! Instantiations must implement functions for PDF and alpha_s access.
   //!

public:
   typedef double(*mu_func)(double,double);

   fastNLOReader(std::string filename);
   fastNLOReader(const fastNLOTable&);
   fastNLOReader(const fastNLOReader&);
   virtual ~fastNLOReader();
   void SetFilename(std::string filename) ;
   void InitScalevariation();
   void SetUnits(fastNLO::EUnits Unit);
   /// Set contribution Id On/Off. Check for Id of a particular contribution with ContrId(...) or use ActivateContribution(...).
   bool SetContributionON(fastNLO::ESMCalculation eCalc , unsigned int Id , bool SetOn = true);
   /// Activate first found contribution of type eCalc and order eOrd
   bool ActivateContribution(fastNLO::ESMCalculation eCalc , fastNLO::ESMOrder eOrd , bool SetOn = true);
   /// Find Id in table of contribution of type eCalc and order eOrd
   int ContrId(const fastNLO::ESMCalculation eCalc, const fastNLO::ESMOrder eOrder) const;
   /// Switch on LO and NLO contributions, deactivate other contributions.
   void SetCoefficientUsageDefault();
   /// Get, if this table is a 'flexible-scale' table or not.
   inline bool GetIsFlexibleScaleTable(fastNLOCoeffAddBase* ctest=NULL) const {
      if ( ctest ) return  ctest->GetNScaleDep() >= 3;
      else if ( B_LO() ) return B_LO()->GetIsFlexibleScale();
      else if ( B_NLO() ) return B_NLO()->GetIsFlexibleScale();
      else if ( B_NNLO() ) return B_NNLO()->GetIsFlexibleScale();
      else return false;
   }
   void SelectProcesses( const std::vector< std::pair<int,int> >& proclist );   //!< tries to select the specified subprocesses for calculation. Prints a warning on failure.
   void SelectProcesses( const std::string& processes, bool symmetric = true );                      //!< tries to select the specified subprocesses for calculation. Prints a warning on failure.

   // ---- setters for specific options ---- //
   void SetNewSqrtS(double NewSqrtS, double OldSqrtS=0 );

   // ---- setters for scales of MuVar tables ---- //
   void SetMuRFunctionalForm(fastNLO::EScaleFunctionalForm func);                       //!< Set the functional form of Mu_R
   void SetMuFFunctionalForm(fastNLO::EScaleFunctionalForm func);                       //!< Set the functional form of Mu_F
   void SetFunctionalForm(fastNLO::EScaleFunctionalForm func , fastNLO::EMuX kMuX);     //!< Set functional form of MuX
   bool SetScaleFactorsMuRMuF(double xmur, double xmuf);                                //!< Set scale factors for MuR and MuF
   void SetExternalFuncForMuR(mu_func);                                                 //!< Set external function for scale calculation (optional)
   void SetExternalFuncForMuF(mu_func);                                                 //!< Set external function for scale calculation (optional)
   void SetExternalConstantForMuR(double MuR);                                          //!< Set value for mu_r if mu_r is chosen to be a constant value (i.e. m_t, or m_Z)
   void SetExternalConstantForMuF(double MuF);                                          //!< Set value for mu_f if mu_f is chosen to be a constant value (i.e. m_t, or m_Z)

   void UseHoppetScaleVariations(bool);

   // ---- Pdf interface ---- //
   void FillPDFCache(double chksum=0., bool lForce=false);                              //!< Prepare for recalculation of cross section with 'new'/updated pdf.
   std::vector<double> GetXFXSqrtS(double x, double muf);                               //!< Interface to GetXFX, but for 'reweighted' sqrt(s)
   virtual std::vector<double> GetXFX(double x, double muf) const = 0;

   // virtual functions for the user interface
   virtual bool InitPDF() = 0;
   virtual double EvolveAlphas(double Q) const = 0;

   // ---- alphas cache ---- //
   void FillAlphasCache(bool lForce=false);                                             //!< prepare for recalculation of cross section with new alpha_s value.

   // --- cache ---- //
   void ResetCache() { fPDFCached=0; fAlphasCached=0;}

   // ---- Do the cross section calculation ---- //
   void CalcCrossSection();
   double RescaleCrossSectionUnits(double binsize, int xunits);                         // Rescale according to kAbsoluteUnits and Ipublunits settings

   std::vector < double > GetCrossSection(bool lNorm = false);      //!< Return vector with all cross section values, normalize on request
   std::vector < double > GetUncertainty(bool lNorm = false);       //!< Return vector with additional uncertainty of cross section values, normalise on request (NOT YET IMPLEMENTED)
   std::vector < double > GetNormCrossSection();                    //!< Return vector with all normalized cross section values
   std::vector < std::map< double, double > > GetCrossSection_vs_x1(); //! Cross section vs. x1 ( XSection_vsX1[bin][<x,xs>] )
   std::vector < std::map< double, double > > GetCrossSection_vs_x2(); //! Cross section vs. x2 ( XSection_vsX1[bin][<x,xs>] )

   std::vector < double > GetReferenceCrossSection();
   std::vector < double > GetQScales();   //!< Order (power of alpha_s) rel. to LO: 0 --> LO, 1 --> NLO
   std::vector < std::vector < double > > GetCrossSection2Dim();


   //! Return struct with vectors containing the cross section values and the selected uncertainty
   XsUncertainty GetScaleUncertainty( const fastNLO::EScaleUncertaintyStyle eScaleUnc, bool lNorm = false);
   XsUncertainty GetAddUncertainty( const fastNLO::EAddUncertaintyStyle eAddUnc, bool lNorm = false);
   //! Function for use with pyext (TODO: Clean this up)
   std::vector< std::vector<double> > GetScaleUncertaintyVec( const fastNLO::EScaleUncertaintyStyle eScaleUnc );
   std::vector< std::vector<double> > GetAddUncertaintyVec( const fastNLO::EAddUncertaintyStyle eAddUnc );

   // ---- Getters for fastNLOReader member variables ---- //
   fastNLO::EScaleFunctionalForm GetMuRFunctionalForm() const { return fMuRFunc; };
   fastNLO::EScaleFunctionalForm GetMuFFunctionalForm() const { return fMuFFunc; };
   fastNLO::EUnits GetUnits() const { return fUnits; };
   mu_func GetExternalFuncForMuR() { return Fct_MuR; };
   mu_func GetExternalFuncForMuF() { return Fct_MuF; };
   double fConst_MuR; //!< Constant _value_ for the renormalization scale. Used only for flexible-scale tables and if requested.
   double fConst_MuF; //!< Constant _value_ for the factorization scale. Used only for flexible-scale tables and if requested.

   double GetScaleFactorMuR() const { return fScaleFacMuR; };
   double GetScaleFactorMuF() const { return fScaleFacMuF; };
   int GetScaleVariation() const { return fScalevar; };
   std::string GetScaleDescription(const fastNLO::ESMOrder eOrder, int iScale=0) const;
   double GetNevt(const fastNLO::ESMOrder eOrder) const;                //!< Get number of events in contribution
   int GetNSubproc(const fastNLO::ESMOrder eOrder) const;               //!< Get number of subprocesses in this contribution
   std::vector < std::vector < std::pair < int,int > > > GetSubprocIndices(const fastNLO::ESMOrder eOrder) const; //!< Get information on the members of each subprocess. Each member of the [iSubproc][iPartonPair] pair is a pair of PDGIds indicating the particles involved in the subprocess.

   int GetNScaleVariations() const;                                                     //!< Get number of available scale variations
   std::vector < double > GetScaleFactors() const;                                           //!< Get list of available scale factors

   // ---- Print outs ---- //
   ///  Print basic info about fastNLO table and its contributions
   void Print(int iprint) const;
   void PrintContributionSummary(int iprint) const;
   void PrintCrossSections() const; //!<  Print cross sections (optimized for double-differential tables)
   void PrintCrossSectionsWithReference();

   // ---- Test virtual functions for reasonable values. ---- //
   bool TestXFX();                                                                      //!< Test if XFX reasonable values
   bool TestAlphas();                                                                   //!< Test if EvolvaAlphas returns a reasonable value


protected:
   fastNLOReader();
   void OrderCoefficients() ;
   //void ReadTable();
   void StripWhitespace(std::string* s);

   void PrintScaleSettings(fastNLO::EMuX kMuX=fastNLO::kMuR);
   void FillBlockBPDFLCsDISv20(fastNLOCoeffAddFix* B);
   void FillBlockBPDFLCsDISv21(fastNLOCoeffAddFlex* B, fastNLOCoeffAddFlex* B0=NULL);
   void FillBlockBPDFLCsHHCv20(fastNLOCoeffAddFix* B);
   void FillBlockBPDFLCsHHCv21(fastNLOCoeffAddFlex* B);
   void CalcAposterioriScaleVariationMuR();
   void CalcAposterioriScaleVariationMuF();
   void FillAlphasCacheInBlockBv20(fastNLOCoeffAddFix* B);
   void FillAlphasCacheInBlockBv21(fastNLOCoeffAddFlex* B);
   double CalcAlphas(double Q);
   double CalcReferenceAlphas();
   double CalcNewPDFChecksum();
   double CalcChecksum(double mu);
   bool PrepareCache();

   void CalcReferenceCrossSection();

   double CalcMu(fastNLO::EMuX kMuX, double scale1 , double scale2 , double scalefactor);
   double FuncMixedOver1(double scale1 , double scale2) ;
   double FuncMixedOver2(double scale1 , double scale2) ;
   double FuncMixedOver4(double scale1 , double scale2) ;
   double FuncMixed2s2Ov2(double scale1 , double scale2) ;
   double FuncMixed2s2Ov4(double scale1 , double scale2) ;
   double FuncPow4Sum(double scale1 , double scale2) ;
   double FuncWgtAvg(double scale1 , double scale2) ;
   double FuncLinearMean(double scale1 , double scale2) ;
   double FuncLinearSum(double scale1 , double scale2) ;
   double FuncMax(double scale1 , double scale2) ;
   double FuncMin(double scale1 , double scale2) ;
   double FuncProd(double scale1 , double scale2) ;
   double FuncExpProd2(double scale1 , double scale2) ;

   void CalcCrossSectionv20(fastNLOCoeffAddFix*  B);
   void CalcCrossSectionv21(fastNLOCoeffAddFlex* B);

   fastNLOCoeffAddBase* B_LO() const {
      //if ( BBlocksSMCalc[fastNLO::kFixedOrder][fastNLO::kLeading] !=0 )
      return (fastNLOCoeffAddBase*) BBlocksSMCalc[fastNLO::kFixedOrder][fastNLO::kLeading];
      // else if ( B_NLO()!= NULL ) return B_NLO();
      // else if ( B_NNLO()!= NULL ) return B_NNLO();
   };
   fastNLOCoeffAddBase* B_NLO() const {
      return (fastNLOCoeffAddBase*) BBlocksSMCalc[fastNLO::kFixedOrder][fastNLO::kNextToLeading];
   };
   fastNLOCoeffAddBase* B_NNLO() const {
      return (fastNLOCoeffAddBase*) BBlocksSMCalc[fastNLO::kFixedOrder][fastNLO::kNextToNextToLeading];
   };
   fastNLOCoeffBase* B_ThC(int n=0) {
      if (BBlocksSMCalc[fastNLO::kThresholdCorrection].empty()) return NULL;
      else return BBlocksSMCalc[fastNLO::kThresholdCorrection][n];
   };
   fastNLOCoeffAddBase* B_Any() const {
      if (B_LO() != NULL ) return B_LO();
      else if ( B_NLO()!= NULL ) return B_NLO();
      else if ( B_NNLO()!= NULL ) return B_NNLO();
      // else if ( B_ThC(0)!= NULL ) return B_ThC(0);
      // else if ( B_ThC(1)!= NULL ) return B_ThC(1);
      else {
         std::cerr<<"Error. Cannot get any additive contribution, but requested."<<std::endl;
         exit(3);
         return NULL;
      }
   };

   // ---- setters for scale variation in v2.0 tables  ---- //
   bool SetScaleVariation(int scalevar);                       //!< Choose the MuF scale variation table

   // ---- human readable strings ---- //
   //static const std::string fContrName[20];
   //static const std::string fOrdName[4][4];
   //static const std::string fNSDep[6];

   bool UpdateProcesses(); //!< Checks if the choosen processes in fselect_processes are compatible to all selected contributions and activate them. Returns true on success false on failure.

protected:
   std::string ffilename;
   int fScalevar;
   double fScaleFacMuR;
   double fScaleFacMuF;
   fastNLO::EScaleFunctionalForm fMuRFunc;
   fastNLO::EScaleFunctionalForm fMuFFunc;
   fastNLO::EUnits               fUnits;
   bool fPDFSuccess;
   double fPDFCached;
   double fAlphasCached;
   mu_func Fct_MuR;                                                                     //!< Function, if you define your functional form for your scale external
   mu_func Fct_MuF;                                                                     //!< Function, if you define your functional form for your scale external

   bool fUseHoppet;
   double fSqrtSovSP = 1;  //!< Center-of-mass 'reweighting'

   std::vector < std::pair<int,int> >* fselected_processes = NULL;   //!< selected processes. When NULL, all processes are used in the calculation

   // ---- pointers to coefftables in fCoeff ---- //
   //    std::vector< std::vector < fastNLOCoeffAddBase* > > fCoAdd;
   //    std::vector< std::vector < fastNLOCoeffMult* > > fCoMult;
   std::vector < std::vector < fastNLOCoeffBase* > > BBlocksSMCalc;                               //!< BlockB's for SM corrections

   // ---- Cross sections ---- //
   std::vector < double > XSection_LO;
   std::vector < double > XSection;
   std::vector < double > dXSection;  // uncertainty sum of x section
   std::vector < double > dX2Section; // squared uncertainty sum of x section
   std::vector < double > QScale_LO;
   std::vector < double > QScale;
   std::vector < std::map< double, double > > fXSection_vsX1; //! Cross section vs. x ( XSection_vsX1[bin][<x,xs>] )
   std::vector < std::map< double, double > > fXSection_vsX2;
   std::vector < std::map< double, double > > fXSection_vsQ2; //RADEK add

   // ----  reference tables ---- //
   std::vector < double > XSectionRef;
   std::vector < double > XSectionRefMixed;
   std::vector < double > XSectionRef_s1;
   std::vector < double > XSectionRef_s2;


};
#endif
