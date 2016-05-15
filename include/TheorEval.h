/**
  @class TheorEval

  @brief Theory expression evaluation 

  This class performs evaluation of expression that describes theoretical
  prediction for a given dataset. The expression may contain APPLgrid terms,
  k-factor terms and digital numbers, all to be read from corresponding sources
  and managed accordingly. After reading, it is converted to reverse polish
  notation and evaluated for every iteration during the fit.

  @author A.Sapronov <sapronov@ifh.de>

  @version 1.2

  @date 2013/08/12
 */
 /*
  2014/03/14: Adding normalisation grids option and more documentation -- AS.
  2015/10/22: Adding fastNLO interface -- DB.
 */
#ifndef TheorEval_h
#define TheorEval_h 1

#include <list>
#include <map>
#include <string>
#include <valarray>
#include <vector>

#ifdef ENABLE_CUTE
#include "CuteInterface.h"
#endif

#ifdef ENABLE_DYTURBO
#include "DyturboInterface.h"
#endif

#include "appl_grid/appl_grid.h"
#include "CommonGrid.h"
//#include "appl_grid/appl_grid.h"

using std::valarray;
using std::vector;
using std::string;
using std::map;

//! Arithmetic token struct
/**
  @struct tToken
  A token to store expression terms. The operator field have precedence
  information:
  0 -- value token, variable or number
  -1, -2  -- brackets (, )
  1 -- operators +, -
  3 -- operators *, /
  4 -- functions sum, avg 
 */
struct tToken {
  short int opr;  // operator flag, includes precedence
  string name;     // string token
  valarray<double> *val; // value token
};

class TheorEval{
 private:
  //! The default constructor is hidden
  TheorEval(){};
 public:
  ~TheorEval();
  //! Proper constructor for TheorEval class
  /*! 
    \param dsID dataset ID for which the expression is evaluated
    \param nTerms the number of terms in the expression
    \param stn array of strings with term names
    \param stt array of strings with term types
    \param sti array of strings with term info
    \param sts array of strings with term sources
    \param ste string with the expression itself

    TheorEval constructor that should be used in the code.
  */
  TheorEval(const int dsID, const int nTerms, const std::vector<string> stn, const std::vector<string> stt,
            const std::vector<string> sti, const std::vector<string> sts, const string& ste);
  
  //! Evaluates array of predictions for requested expression
  /*!
    \param iorder order for grids evaluation
    \param mur renormalisation scale used in grid evaluation
    \param muf factorisation scale
    \param vte the resulting valarray with predictions

    The prupose of this method is to evaluate the expression for current
    iteration. It updates the expression components and folds the reverse
    polish notation to calculate the result.
   */
  int Evaluate(valarray<double> &vte );

  //! Set custom CKM matrix for APPLgrid
  /*!
    \param v_ckm the 3x3 CMK matrix to be set in the APPLgrid

    The CKM matrix values used in APPLgrids calculations can be updated
    here.
   */
  int setCKM(const vector<double> &v_ckm);

  //! Set dataset bins
  /*!
    \param nBinDim the binning dimension (only 1d is supported at the moment)
    \param nPoints number of points (bins)
    \param binFlags array with flags for each bin
    \param allBins array of bin boundaries

    This method sets the binning of the dataset.
   */
  int setBins(int nBinDim, int nPoints, int *binFlags, double *allBins);
  //! Initializes sources for theoretical predictions
  /*!
   After the datasets with expressions are read, this method initialises
   terms such as applgrids and k-factor tabels from their sources.
   */
  int initTheory();
  //! Returns numebr of bins in the current dataset
  int getNbins();
  //! Returs vector with bin flags for current dataset
  const vector<int> *getBinFlags() const { return &_binFlags; }
  //! Selects if we have a proton-antiproton collision
  void SetCollisions(int ppbar) {_ppbar = (ppbar == 1);};
  void SetDynamicScale(float dynscale) {_dynamicscale = dynscale;};
  void SetYLow(float y) {_ylow = y;};
  void SetYHigh(float y) {_yhigh = y;};
  void SetMLow(float m) {_mlow = m;};
  void SetMHigh(float m) {_mhigh = m;};
  void SetNormalised(int normalised) {_normalised = (normalised == 1);};
  void SetBinDensity(int bindensity) {_bindensity = bindensity;};
  void SetMurMufDef(int MurDef, int MufDef) { _MurDef = MurDef; _MufDef = MufDef;}; //!< Set mur and muf definition for fastNLO flexible-scale tables
  void SetOrdScales(int iord, double mur, double muf, double mures, double muC3 = 1.) { _iOrd=iord; _xmur=mur; _xmuf=muf; _xmures = mures; _xmuC3 = muC3;}; //!< set order and scale factors
  void GetOrdScales(int &iord, double &mur, double &muf, double &mures) { iord=_iOrd; mur=_xmur; muf=_xmuf; mures = _xmures;}; //!< get order and scale factors
  void GetOrdScales(int &iord, double &mur, double &muf, double &mures, double &muC3) { iord=_iOrd; mur=_xmur; muf=_xmuf; mures = _xmures; muC3 = _xmuC3;}; //!< get order and scale factors
  void ChangeTheorySource(string term, string source);
  string GetTheorySource(string term);

 private:
  //! Checks that the bin boundaries in theory sources are complied with data ones.
  int checkBins();
  //! Tokenize symbolic expression to a list of string tokens
  int tokenize(string &, list<string> &);
  //! The tokens are assigned to corresponding values
  int assignTokens(list<tToken> &);
  //! The expression of tokens is converted to RPN
  int convertToRPN(list<tToken> &);
  //! Initialise terms
  /*!
   Depending on the term type, the corresponding initialization method is called
  */
  int initTerm(int, valarray<double> *);

  //! Initialise applgrid-based term
  /*!
   \param iterm expression term index
   \param val valarray pointer which should be associated with the term

   Initializes the applgrid-based grids and associates the term valarrays with them. 
  */
  int initGridTerm(int iterm, valarray<double> *val);
  //! Initialise K-factor term
  int initKfTerm(int, valarray<double> *);
  void initCuteTerm(int iterm, valarray<double> *val);
  void initDyturboTerm(int iterm, valarray<double> *val);
  //! Get current grid values into the tokens
  int getGridValues();
  void getCuteValues();
  void getDyturboValues();

 private:
  int _dsId;
  int _iOrd;
  double _xmur;
  double _xmuf;
  double _xmures;
  double _xmuC3;
  int _nTerms;
  double _units;
  vector<string> _termNames;
  vector<string> _termTypes;
  vector<string> _termInfos;
  vector<string> _termSources;
  string _expr;
  vector<int> _binFlags;
  vector<vector<double> > _dsBins;
  vector<valarray<double> > _kFactors;

  //CUTE
  //map of cute predictions, is iterated on all terms when evaluating predictions
#ifdef ENABLE_CUTE
  map <Cute*, valarray<double>* > _mapCuteTerms;
#endif

#ifdef ENABLE_DYTURBO
  map <Dyturbo*, valarray<double>* > _mapDyturboTerms;
#endif

  //cute
  //_exprRPN is the vector of all the terms, in tToken format

  /// Reverse polish notation of the expression
  vector<tToken> _exprRPN;
  map<CommonGrid*, valarray<double>* > _mapGridToken;
  map<string, valarray<double>* > _mapInitdTerms;

  /// Normalised theory
  bool _normalised;

  /// ppbar PDF
  bool _ppbar;

  /// fastNLO flexible-scale defintions
  int _MurDef;
  int _MufDef;

  /// bin-by-bin dynamic scale
  float _dynamicscale;
  float _ylow;
  float _yhigh;
  float _mlow;
  float _mhigh;

  //number of calculation points per bin
  int _bindensity;
};

typedef map <int, TheorEval* > tTEmap;

/// global dataset to theory evaluation pointer map
extern tTEmap gTEmap;

#endif
