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


class ReactionTheory;
class TermData;

using std::valarray;
using std::vector;
using std::string;
using std::map;
using std::list;

//TheorEval live in map<int,TheorEval*>gTEmap: datasetID->TheorEval
//TheorEval are created in set_theor_eval_(dataset id), in src/ftheor_eval.cc

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
  5 -- Matrix . Vector / vector . Matrix
 */
struct tToken {
  short int opr;  // operator flag, includes precedence
  string name;     // string token
  valarray<double> *val; // value token
  int narg;
  bool ownsVal=true;//true if val should be destroyed when this tToken is destroyed
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
  void Evaluate(valarray<double> &vte );

  //! Set dataset bins
  /*!
    \param nBinDim the binning dimension (only 1d is supported at the moment)
    \param nPoints number of points (bins)
    \param binFlags array with flags for each bin
    \param allBins array of bin boundaries (nope --Ivan)

    This method sets the binning of the dataset.
   */
  void setBins(int nBinDim, int nPoints, int *binFlags, double *allBins,map<string,size_t>&columnNameMap);
  //! Initializes sources for theoretical predictions
  /*!
   After the datasets with expressions are read, this method initialises
   terms such as applgrids and k-factor tabels from their sources.
   */
  void initTheory();
  //The following methods provide access to bins of this dataset
  size_t getNbins()const{return _dsBins[0].size();}//return number of bins (datapoints)
  bool hasBinColumn(const string&columnName)const{return columnNameMap.count(columnName)==1;}
  const valarray<double>* getBinColumn(const string& columnName)const;//return nullptr if not found
  const vector<int>* getBinFlags()const{return&_binFlags;}
  void SetNormalised(int normalised) {_normalised = (normalised == 1);};
  void SetNormalisation(double normalisation) {_normalisation = normalisation;};

  //The following methods are used by chi2scan and allow changing a theory input file
  void ChangeTheorySource(string term, string source);
  string GetTheorySource(string term);

 private:
  //! Checks that the bin boundaries in theory sources are complied with data ones.
  int checkBins();
  //! Tokenize symbolic expression to a list of string tokens
  int tokenize(string &, list<string> &);
  //! The tokens are assigned to corresponding values
  void assignTokens(list<tToken> &);
  //! The expression of tokens is converted to RPN
  void convertToRPN(list<tToken> &);
  //! Initialise terms
  /*!
   Depending on the term type, the corresponding initialization method is called
  */
  void initTerm(int, valarray<double> *);

  /*!
   \brief Initializes a token corresponding to a reaction term
   \param[in,out] token The token to be initialized
   \param[in]     name  The name of the reaction term. Will be assigned as the name of the token.
   \details If the reaction with the given \p name does not exist, issues a fatal error
  */
  void initReactionToken(tToken& token,const string& name);
  //! Initialise reaction term
  void initReactionTerm(int iterm, valarray<double>* val, bool change_source = false);
  //! Update the reaction values into the tokens
  void updateReactionValues();

  int _nTerms;
  vector<string> _termNames;
  vector<string> _termTypes;//we do not need this anymore, the one and only term type is "reaction"
  vector<string> _termInfos;
  vector<string> _termSources;
  string _expr;
  vector<int> _binFlags;//_binFlags[i] is flag for bin i. Flag=0 means bin is disabled and excluded from the fit. Flag=1 means enabled.
  vector<valarray<double> >_dsBins;//_dsBins[i][j] is value in 'Bin'-type column i, row j, as provided in datafile
  map<string,size_t>columnNameMap;//Map from column name to column number i; _dsBins[i] is corresp. column

  /// Reverse polish notation of the expression
  vector<tToken> _exprRPN;
  map<string, valarray<double>* > _mapInitdTerms;

  /// Normalised theory
  bool _normalised=false;
  double _normalisation=1;
  vector<TermData*>term_datas;
public:
  /// also keep some dataset information:
  string _ds_name;///Name
  //TODO: Field names could be better --Ivan
  int _dsId;   //Id of dataset, which is the number of corresponding entry in steering.txt
  int _dsIndex;//Index of dataset, as given in the datafile
};
typedef map <int, TheorEval* > tTEmap;
typedef map <string, string> tReactionLibsmap;
typedef map <string, ReactionTheory *> tNameReactionmap;
/// global dataset to theory evaluation pointer map
extern tTEmap gTEmap;

/// global map of reaction libraries
extern tReactionLibsmap gReactionLibs;

/// global map of reaction names
extern tNameReactionmap gNameReaction;
#endif
