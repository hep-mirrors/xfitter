#ifndef TheorEval_h
#define TheorEval_h 1

#include <list>
#include <map>
#include <string>
#include <valarray>
#include <vector>

#include "appl_grid/appl_grid.h"

using namespace std;

struct tToken {
  short int opr;  // operator flag, includes precedence
  string name;     // string token
  valarray<double> *val; // value token
};

//! Manages theretical evaluations with APPLgrid, FastNLO, k-factors
/*!
 * The class instance correspond to one dataset and performs operations
 * related to the grid managment: open/close, check bins, calculate
 * cross sections, apply cuts, etc.
 */
class TheorEval{
 private:
  //! The default constructor is hidden
  TheorEval(){};
 public:
  ~TheorEval();
  //! Proper constructor for TheorEval class
  TheorEval(const int dsID, const int nTerms, const string* stn, const string* stt, 
            const string* sts, const string& ste);
  
  //! Evaluates array of predictions for requested expression
  int Evaluate(const int iorder, const double mur, const double muf, valarray<double> &vte );
  //! Set custom CKM matrix for APPLgrid
  int setCKM(const vector<double> &);
  //! Set dataset bins
  int setBins(int nBinDim, int nPoints, int *binFlags, double *allBins);
  //! Initializes sources for theoretical predictions
  int initTheory();
  int getNbins();
  int setUnits(double units){ _units = units;};
  const vector<int> *getBinFlags() const { return &_binFlags; }
  void SetCollisions(int ppbar) {_ppbar = (ppbar == 1);};
  void SetNormalised(int normalised) {_normalised = (normalised == 1);};
  void SetDynamicScale(float dynscale) {_dynamicscale = dynscale;};

 private:
  //! Checks that the bin boundaries are complied with data ones.
  int checkBins();
  //! Tokenize symbolic expression to a list of string tokens
  int tokenize(string &, list<string> &);
  int assignTokens(list<tToken> &);
  int convertToRPN(list<tToken> &);
  int initTerm(int, valarray<double> *);
  int initGridTerm(int, valarray<double> *);
  int initKfTerm(int, valarray<double> *);
  //! Get current grid values into the tokens
  int getGridValues(const int iorder, const double mur, const double muf);

 private:
  int _dsId;
  int _nTerms;
  double _units;
  vector<string> _termNames;
  vector<string> _termTypes;
  vector<string> _termSources;
  string _expr;
  vector<int> _binFlags;
  vector<vector<double> > _dsBins;
  vector<valarray<double> > _kFactors;

  /// Reverse polish notation of the expression
  vector<tToken> _exprRPN;
  map<appl::grid*, valarray<double>* > _mapGridToken;
  map<string, valarray<double>* > _mapInitdTerms;

  //ppbar PDF
  bool _ppbar;
  //Normalised cross section
  bool _normalised;
  //bin-by-bin dynamic scale
  float _dynamicscale;
};

#endif
