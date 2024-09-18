#pragma once
#include<cstddef>
using std::size_t;
#include<string>
using std::string;
#include<vector>
using std::vector;
#include<map>
using std::map;
#include<valarray>
using std::valarray;
#include"TheorEval.h"
class ReactionTheory;
namespace xfitter{class BaseEvolution;}
/*
TermData is a class that provides interface to all parameters for a given reaction term
TermData is responsible for:
*Providing interface to reaction parameters, resolving overwrites and global vs specific parameters
*Reading special parameters "evolution", "evolution1" and "evolution2", get pointer[s] to correct evolution[s]
*Switch XFXlike and other wrappers

One instance of TermData exists for each reaction term.
Instances of TheorEval manage all instances of TermData and pass them to ReactionTheory

It is intended to be used by ReactionTheory
*/
class TermData{
public:
  TermData(unsigned id,ReactionTheory*,TheorEval*parent,const char*TermInfo);
  //Unique id of this term
  //In the past this was calculated as 1000*datasetID+number_of_term (see TheorEval::initReactionTerm)
  //In ReactionTheory called it was also incorrectly called dataSetID
  const unsigned id;
  //Return true if parameter with given name exists, false otherwise
  bool         hasParam (const string&parameterName);
  //The following 3 methods return reaction parameter by its name, taking into account that term-specific parameters overshadow global etc.
  //They either return a parameter of the requested type, or issue a fatal error
  const double*getParamD(const string&parameterName);//returns pointer because double parameters can change at each iteration
  int          getParamI(const string&parameterName);
  string       getParamS(const string&parameterName);
  ReactionTheory*reaction;
  void actualizeWrappers();//see wrappers below
  xfitter::BaseEvolution*getPDF(int i=0);//i is either 0 for evolution1, or 1 for evolution2
  //The following methods are used to get binning information for the dataset of this term
  size_t getNbins(){return parent->getNbins();}
  bool                   hasBinColumn(const string&columnName){return parent->hasBinColumn(columnName);}
  const valarray<double>&getBinColumn(const string&columnName);
  const valarray<double>*getBinColumnOrNull(const string&columnName){return parent->getBinColumn(columnName);}//same, but return nullptr instead of issuing an error if the column does not exist
  const vector<int>&     getBinFlags(){return*parent->getBinFlags();}//0 means bin is disabled, 1 means enabled. Disabled bins are excluded from the fit
  //The following pointer can be used by ReactionTheory to store some additional data
  //for each reaction term. It should be managed by ReactioTheory only, do not touch it from elsewhere
  void*reactionData=nullptr;
  //array that the reaction should fill with predictions at compute()
  //This is intended to be used by TheorEval
  //val might be moved somewhre else in the future, please do not use it outside of TheorEval
  //this is a "weak" pointer, the actual valarray is owned by TheorEval
  valarray<double>*val=nullptr;
  int _ncpu; // number of parallel threads
private:
  TheorEval*parent;//The instance of TheorEval that manages this instance of TermData
  map<string,string>term_info;//Map key->value
  int getStringFromTermOrReaction(const string&,string*);
};
/* Wrappers
  If a ReactionTheory uses these wrapper, it must call TermData::actualizeWrappers
  for each term at each iteration, to make sure the wrappers wrap the correct evolution
*/
extern "C"{
void pdf_xfxq_wrapper_ (const double& x, const double& Q,double* results);
void pdf_xfxq_wrapper1_(const double& x, const double& Q,double* results);
double pdf_ixfxq_wrapper_ (const int& id, const double& x, const double& Q);
double pdf_ixfxq_wrapper1_(const int& id, const double& x, const double& Q);

double alphas_wrapper_(const double& Q);
}
