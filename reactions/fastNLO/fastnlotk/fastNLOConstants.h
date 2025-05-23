#ifndef __fnloconstants__
#define __fnloconstants__

// NEVER EVER include a project's internal config.h in installable header files!
// Use for conditional compilation only in .cc source code files.
// Otherwise conflicts with other linked projects are to be expected.
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef FNLO_NAME
#define FNLO_NAME       "fastNLO_toolkit"
#define FNLO_SUBPROJECT "toolkit"
#define FNLO_VERSION    "2.6.0"
#define FNLO_GITREV     "3067-3-g1f8468f0"
#define FNLO_AUTHORS    "D. Britzger, T. Kluge, K. Rabbertz, F. Stober, G. Sieber, M. Wobisch"
#define FNLO_WEBPAGE    "http://projects.hepforge.org/fastnlo"
#define FNLO_AUTHORSv14 "T. Kluge, K. Rabbertz, M. Wobisch"
#define FNLO_QUOTEv14   "Proc. DIS 2006, 483 (2006), hep-ph/0609285."
#define FNLO_AUTHORSv2  "D. Britzger, K. Rabbertz, F. Stober, M. Wobisch"
#define FNLO_QUOTEv2    "Proc. DIS 2012, 217 (2012), arXiv:1208.3641."
#define FNLO_YEARS      "2005-2025"
#endif

// KR: Replace by precompiler defines
const double TWOPI    =  6.28318530717958647692528;
const double TWOPISQR = 39.47841760435743447533796;
const double TOCL90   =  1.64485362695147271486385; // SQRT(2.D0)*InvERF(0.9D0)
//#define TWOPI (2.*M_PI)
//#define TWOPISQR (4.*M_PI*M_PI)
// PDG values 2017, MSbar except MT
#define PDG_MU   (0.0022)
#define PDG_MD   (0.0047)
#define PDG_MS   (0.096)
#define PDG_MC   (1.28)
#define PDG_MB   (4.18)
#define PDG_MT   (173.1)
#define PDG_MZ   (91.1876)
#define PDG_ASMZ (0.1182)

namespace fastNLO {

   // ---- typedefs ---- //
   typedef std::vector<double > v1d;
   typedef std::vector<std::vector<double > > v2d;
   typedef std::vector<std::vector<std::vector<double > > > v3d;
   typedef std::vector<std::vector<std::vector<std::vector<double > > > > v4d;
   typedef std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > > v5d;
   typedef std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > > > v6d;
   typedef std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > > > > v7d;

   // ---- constants ---- //
   static const std::set<int> CompatibleVersions{20000,21000,22000,23000,23500,23600,25000,26000};
   const int tabversion   = 26000;
   const int tablemagicno = 1234567890;
   // separating character between entries in table
   const char sep[] = "\n";
   //   const double TWOPI = 6.28318530717958647692528;


   // ---- enumerators ---- //
   enum EMuX {
      kMuR                        = 0,    // renormalization scale
      kMuF                        = 1     // factorization scale
   };

   enum EScaleFunctionalForm {
      kScale1                   = 0,    // e.g. mu^2 = Q^2
      kScale2                   = 1,    // e.g. mu^2 = pt^2
      kQuadraticSum             = 2,    // e.g. mu^2 = ( Q^2 + pt^2 )
      kQuadraticMean            = 3,    // e.g. mu^2 = ( Q^2 + pt^2 ) / 2
      kQuadraticSumOver4        = 4,    // e.g. mu^2 = ( Q^2 + pt^2 ) / 4
      kLinearMean               = 5,    // e.g. mu^2 = (( Q + pt ) / 2 )^2
      kLinearSum                = 6,    // e.g. mu^2 = (( Q + pt ))^2
      kScaleMax                 = 7,    // e.g. mu^2 = max( Q^2, pt^2)
      kScaleMin                 = 8,    // e.g. mu^2 = min( Q^2, pt^2)
      kProd                     = 9,    // e.g. mu^2 = (scale1 * scale2) ^2
      kS2plusS1half             = 10,   // e.g. mu^2 = (scale1^1/2 + scale2^2)
      kPow4Sum                  = 11,   // e.g. mu^2 = sqrt((scale1^4 + scale2^4))
      kWgtAvg                   = 12,   // e.g. mu^2 = sqrt( (scale1^4 + scale2^4)/(scale1^2 + scale2^2)) [weighted average]
      kS2plusS1fourth           = 13,   // e.g. mu^2 = (scale1^1/4 + scale2^2)
      kExpProd2                 = 14,   // e.g. mu^2 = (scale1 * exp(0.3 * scale2)) ^2
      kExtern                   = 15,   // define an external function for your scale
      kConst                    = 16,   // e.g. mu^2 = c, while c is a constant and could be for instance the top-mass
   };

   enum ESMCalculation {
      kFixedOrder               = 0,    // Fixed order calculation (pQCD)
      kThresholdCorrection      = 1,    // Threshold corrections
      kElectroWeakCorrection    = 2,    // Electroweak corrections
      kNonPerturbativeCorrection= 3,    // Non-perturbative corrections|Hadronisation corrections
      kContactInteraction       = 10    // Contact interactions
   };

   enum ESMOrder {
      kLeading                  = 0,    // LO,   1-loop, LO MC
      kNextToLeading            = 1,    // NLO,  2-loop, NLO MC
      kNextToNextToLeading      = 2     // NNLO, 3-loop, NNLO MC
   };

   enum EUnits {
      kAbsoluteUnits            = 0,    // calculate the cross section in barn for each publicated bin
      kPublicationUnits         = 1     // calculate the cross section in units as given in the according publication
   };

   enum EScaleUncertaintyStyle {
      kScaleNone                = 0,    // no scale uncertainty, only central scale (mu_r,mu_f) = (1,1) evaluated
      kSymmetricTwoPoint        = 1,    // symmetric (mu_r,mu_f) scale variations by factors (1/2,1/2), (2,2)
      kAsymmetricSixPoint       = 2,    // asymmetric (mu_r,mu_f) scale variations by factors (1/2,1/2), (2,2) plus
                                        // (1/2,1), (1,1/2), (1,2), (2,1)
      kAsymmetricRatio          = 3     // 30(+1) combinations for mu_r, mu_f in numerator and denominator
   };

   constexpr std::string_view ScaleUncertaintyStyle_to_string(EScaleUncertaintyStyle estyle) {
      switch (estyle) {
         case kScaleNone:          return "NN";
         case kSymmetricTwoPoint:  return "2P";
         case kAsymmetricSixPoint: return "6P";
         case kAsymmetricRatio:    return "30";
         default:                  return "??";
      }
   }

   constexpr std::optional<EScaleUncertaintyStyle> ScaleUncertaintyStyle_to_enum(std::string_view sstyle) {
      if (sstyle == "NN") return kScaleNone;
      if (sstyle == "2P") return kSymmetricTwoPoint;
      if (sstyle == "6P") return kAsymmetricSixPoint;
      if (sstyle == "30") return kAsymmetricRatio;
      return{kScaleNone};
   }

   enum ENumUncertaintyStyle {
      kNumNone                  = 0,    // no numerical uncertainty
      kStatInt                  = 1,    // statistical uncertainty from MC integration
      kApproxBias               = 2     // bias from approximation through interpolation
   };

   constexpr std::string_view NumUncertaintyStyle_to_string(ENumUncertaintyStyle estyle) {
      switch (estyle) {
         case kNumNone:    return "NN";
         case kStatInt:    return "ST";
         case kApproxBias: return "AB";
         default:          return "??";
      }
   }

   constexpr std::optional<ENumUncertaintyStyle> NumUncertaintyStyle_to_enum(std::string_view sstyle) {
      if (sstyle == "NN") return kNumNone;
      if (sstyle == "ST") return kStatInt;
      if (sstyle == "AB") return kApproxBias;
      return{kNumNone};
   }

   enum EPDFUncertaintyStyle {
      kPDFNone                  = 0,    // No PDF uncertainty, only averaged cross section result evaluated (Correct for NNPDF, wrong otherwise!)
      kLHAPDF6                  = 1,    // LHAPDF6 uncertainties (recommended if LHAPDF6 is available)
      kHessianSymmetric         = 2,    // symmetric Hessian PDF uncertainties (ABM)
      kHessianAsymmetric        = 3,    // asymmetric Hessian PDF uncertainties
      kHessianAsymmetricMax     = 4,    // asymmetric Hessian PDF uncertainties with pairwise max deviations per eigenvector (CTEQ,MRST|MSTW)
      kHessianCTEQCL68          = 5,    // like kHessianAsymmetricMax, but with uncertainties rescaled to CL68
      kMCSampling               = 6,    // statistical sampling PDF uncertainties (and central value) (NNPDF)
      kHeraPDF10                = 7     // HERAPDF 1.0 uncertainties
   };

   constexpr std::string_view PDFUncertaintyStyle_to_string(EPDFUncertaintyStyle estyle) {
      switch (estyle) {
         case kPDFNone:              return "NN";
         case kLHAPDF6:              return "L6";
         case kHessianSymmetric:     return "HS";
         case kHessianAsymmetric:    return "HA";
         case kHessianAsymmetricMax: return "HP";
         case kHessianCTEQCL68:      return "HC";
         case kMCSampling:           return "MC";
         case kHeraPDF10:            return "H1"; // not yet implemented
         default:                    return "??";
      }
   }

   constexpr std::optional<EPDFUncertaintyStyle> PDFUncertaintyStyle_to_enum(std::string_view sstyle) {
      if (sstyle == "NN") return kPDFNone;
      if (sstyle == "L6") return kLHAPDF6;
      if (sstyle == "HS") return kHessianSymmetric;
      if (sstyle == "HA") return kHessianAsymmetric;
      if (sstyle == "HP") return kHessianAsymmetricMax;
      if (sstyle == "HC") return kHessianCTEQCL68;
      if (sstyle == "MC") return kMCSampling;
      if (sstyle == "H1") return kHeraPDF10;
      return{kPDFNone};
   }

   enum EAsUncertaintyStyle {
      kAsNone                   = 0,    // no a_s uncertainty
      kAsGRV                    = 1,    // a_s(M_Z) uncertainty with GRV evolution
   };

   constexpr std::string_view AsUncertaintyStyle_to_string(EAsUncertaintyStyle estyle) {
      switch (estyle) {
         case kAsNone:             return "NN";
         case kAsGRV:              return "AS";
         default:                  return "??";
      }
   }

   constexpr std::optional<EAsUncertaintyStyle> AsUncertaintyStyle_to_enum(std::string_view sstyle) {
      if (sstyle == "NN") return kAsNone;
      if (sstyle == "AS") return kAsGRV;
      return{kAsNone};
   }

   enum EMerge {  //!< mergeing options.
      kMerge, //!< Calculate weighted average (default. Nevt usually set externally by generator code).
      kAdd, //!< Add (Append)! Do not merge, but add two tables together (fully unweighted) (1+1=2).
      kUnweighted, //!< Calculated unweighted average (usually better: take kNumEvent).
      kAttach, //!< Add (Append)! Same functionality as 'add' but subprocesses are attached and file size increases.
      kNumEvent, kNumEventBinProc, //!< Calculate weighted average, using w = num entries
      kSumW2,    kSumW2BinProc, //!< Calculate weighted average , using w = sum(weight**2)
      kSumSig2,  kSumSig2BinProc, //!< Calculate weighted average, using w = sum(sig**2) [sig ~ wgt*as*pdf]
      kSumUser,  kSumUserBinProc, //!< Calculate weighted average, using w = sum(sig) [sig ~ wgt*as*pdf], or 'user-specified' weights
      kMedian, kMean, //!< build median or median value of many tables (option not applicable to member function, because many tables are needed as input).
      kUndefined //!< Error
   };


   // ---- some names for nice output ---- //
   const std::string _ContrName[20] = {
      "Fixed order calculation", "Threshold corrections", "Electroweak corrections", "Non-perturbative corrections",
      "Undefined", "Undefined", "Undefined", "Undefined", "Undefined", "Undefined", "Undefined",
      "Quark compositeness", "ADD-LED", "TeV 1-ED", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown"
   };
   const std::string _OrdName[4][4] = {
      { "LO",     "NLO",    "NNLO"   , "N3LO"    },
      { "1-loop", "2-loop", "3-loop" , "4-loop"  },
      { "Undef" , "Undef" , "Undef"  , "Undef"   },
      { "LO MC" , "NLO MC", "NNLO MC", "N3LO MC" }
   };
   const std::string _fNSDep[6] = {"v2.0","v2.0","v2.0","v2.1","v2.2","v2.2"};



   // ---- Some shapes for nice output ---- //
   //
#ifndef SWIG
   const std::string _CSEP20("####################");
   const std::string _DSEP20("====================");
   const std::string _SSEP20("--------------------");
   const std::string _TSEP20(" - - - - - - - - - -");
   const std::string _CSEP20C(" ######################");
   const std::string _DSEP20C(" #=====================");
   const std::string _SSEP20C(" #---------------------");
   const std::string _TSEP20C(" #- - - - - - - - - - -");
   const std::string _CSEP40  = _CSEP20  + _CSEP20;
   const std::string _DSEP40  = _DSEP20  + _DSEP20;
   const std::string _SSEP40  = _SSEP20  + _SSEP20;
   const std::string _TSEP40  = _TSEP20  + _TSEP20;
   const std::string _CSEP40C = _CSEP20C + _CSEP20;
   const std::string _DSEP40C = _DSEP20C + _DSEP20;
   const std::string _SSEP40C = _SSEP20C + _SSEP20;
   const std::string _TSEP40C = _TSEP20C + _TSEP20;
   const std::string _CSEPS  = _CSEP40  + _CSEP40;
   const std::string _DSEPS  = _DSEP40  + _DSEP40;
   const std::string _SSEPS  = _SSEP40  + _SSEP40;
   const std::string _TSEPS  = _TSEP40  + _TSEP40;
   const std::string _CSEPSC = _CSEP40C + _CSEP40;
   const std::string _DSEPSC = _DSEP40C + _DSEP40;
   const std::string _SSEPSC = _SSEP40C + _SSEP40;
   const std::string _TSEPSC = _TSEP40C + _TSEP40;
   const std::string _CSEPL  = _CSEPS  + _CSEPS ;
   const std::string _DSEPL  = _DSEPS  + _DSEPS ;
   const std::string _SSEPL  = _SSEPS  + _SSEPS ;
   const std::string _TSEPL  = _TSEPS  + _TSEPS ;
   const std::string _CSEPLC = _CSEPSC + _CSEPS ;
   const std::string _DSEPLC = _DSEPSC + _DSEPS ;
   const std::string _SSEPLC = _SSEPSC + _SSEPS ;
   const std::string _TSEPLC = _TSEPSC + _TSEPS ;
#endif
}

#endif
