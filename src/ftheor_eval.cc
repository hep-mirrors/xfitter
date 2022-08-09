/*!
 @file theor_eval.cc
 @date Tue Aug 6 2013
 @author Andrey Sapronov <Andrey.Sapronov@cern.ch>

 Contains Fortran interface functions to operate with theoretical
 predictions obtained via fast cross section evaluation methods,
 e.g. APPLgrid,  FastNLO and k-Factors.
 */

#include <vector>
#include <fstream>
#include <valarray>

// #include "get_pdfs.h"
#include "xfitter_cpp.h"
#include "xfitter_cpp_base.h"

#include "TheorEval.h"
#include <yaml-cpp/yaml.h>
#include "ReactionTheory.h"
#include "xfitter_pars.h"
#include"dependent_pars.h"

#include "BasePdfParam.h"
#include "BasePdfDecomposition.h"
#include "BaseEvolution.h"
#include "BaseMinimizer.h"

using namespace std;

extern "C" {
  // ! check consistency of C and Fortran common blocks
  void common_check_(int *i);
}

void common_check_(int *i) {
  if ( *i != steering_.steering_check) {
    string text = "F: Inconsistency of the fortran common steering and C-structure steering_. Check steering.inc and xfitter_cpp.h that the list of variables matches";
    hf_errlog_(17032505,text.c_str(),text.size());
  }
}

extern "C" {
  int set_theor_eval_(int *dsId);//, int *nTerms, char **TermName, char **TermType,
//    char **TermSource, char *TermExpr);
  int set_theor_bins_(int *dsId, int *nBinDimension, int *nPoints, int *binFlags,
                      double *allBins, char binNames[10][80]);
//  int set_theor_units_(int *dsId, double *units);
  void init_theor_eval_(int *dsId);
  void update_theor_ckm_();
  void get_theor_eval_(int *dsId, int* np, int* idx);
  void close_theor_eval_();
  void init_at_iteration_(); ///< Loop over reactions, initialize them
  void fcn3action_();      ///< Loop over reactions, call actionAtFCN3
  void error_band_action_(const int& i); ///< Loop over rections, call error_band_action
}

/// global dataset to theory evaluation pointer map
tTEmap gTEmap;
tNameReactionmap gNameReaction;

const size_t NTERMMAX      =128;
const size_t TERMNAME_LEN  =32;
const size_t TERMTYPE_LEN  =80;
const size_t TERMINFO_LEN  =4096;
const size_t TERMSOURCE_LEN=256;
const size_t THEOREXPR_LEN =10000;
extern struct thexpr_cb {
  double dynscale;
  int nterms;
  char termname  [NTERMMAX][TERMNAME_LEN];
  char termtype  [NTERMMAX][TERMTYPE_LEN];
  char terminfo  [NTERMMAX][TERMINFO_LEN];
  char termsource[NTERMMAX][TERMSOURCE_LEN];
  char theorexpr[THEOREXPR_LEN];
  int normalised;
  double normalisation;
  double datainfo[100];
  int ninfo;  // dataset info as well
  char CInfo[100][80];
  char dsname[80];
  int  ds_index;
} theorexpr_;

extern struct ord_scales {
   double datasetmur[150];
   double datasetmuf[150];
   int datasetiorder[150];
} cscales_;

inline std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}


/*!
 Creates theory evaluation object and adds it to the global map by
 dataset ID.
 write details on argumets
 */
int set_theor_eval_(int *dsId){
  // convert fortran strings to c++
  vector<string> stn(theorexpr_.nterms);
  vector<string> stt(theorexpr_.nterms);
  vector<string> sti(theorexpr_.nterms);
  vector<string> sts(theorexpr_.nterms);
  for ( int i = 0; i< theorexpr_.nterms; i++){
    stn[i]=stringFromFortran(theorexpr_.termname  [i],TERMNAME_LEN);
    stt[i]=stringFromFortran(theorexpr_.termtype  [i],TERMTYPE_LEN);
    sti[i]=stringFromFortran(theorexpr_.terminfo  [i],TERMINFO_LEN);
    sts[i]=stringFromFortran(theorexpr_.termsource[i],TERMSOURCE_LEN);
  }
  string ste=stringFromFortran(theorexpr_.theorexpr,THEOREXPR_LEN);
  TheorEval *te = new TheorEval(*dsId, theorexpr_.nterms, stn, stt, sti, sts, ste);

  /* sometime in xFitter 2.2 CINFO support was dropped
  // Store CINFO
  for (int i=0; i<theorexpr_.ninfo; i++) {
    std::string n(theorexpr_.CInfo[i]);
    n = n.substr(0,80);
    n.erase(std::remove(n.begin(), n.end(), ' '), n.end());
    te->AddDSParameter(n, theorexpr_.datainfo[i]);
  }
  */
  // Store some dataset information
  te->_ds_name=stringFromFortran(theorexpr_.dsname,80);
  te->_dsId=*dsId;
  te->_dsIndex=theorexpr_.ds_index;
  te->SetNormalised(theorexpr_.normalised);
  te->SetNormalisation(theorexpr_.normalisation);

  tTEmap::iterator it = gTEmap.find(*dsId);
  if (it == gTEmap.end() ) { gTEmap[*dsId] = te; }
  else {
    cerr<<"[ERROR] Theory evaluation for dataset ID "<<*dsId<<" already exists."<<endl;
    hf_errlog(19042010,"F: Programming error: TheorEval already exists; see stderr");
  }

  return 1;
}

/*!
 Pass bin information from fortran to instances of TheorEval
 dsId            - dataset ID. Identifies instance of TheorEval
 nBinDimension   - number of bin columns
 nPoints         - number of points in the dataset (=number of rows)
 binFlags[i]     - flag of point i
   Flag=1 means the point is enabled
   Flag=0 means the point is disabled and is excluded from the fit
 allBins[10*j+i] - value at row (datapoint) j in column i (10 is max value of nBinDimension)
 binNames[i]     - name of bin column i
 */
const size_t COLUMN_NAME_LEN=80;
int set_theor_bins_(int *dsId, int *nBinDimension, int *nPoints, int *binFlags,
                    double *allBins, char binNames[10][COLUMN_NAME_LEN])
{
  tTEmap::iterator it = gTEmap.find(*dsId);
  if (it == gTEmap.end() ) {
    cout << "ERROR: Theory evaluation for dataset ID " << *dsId
    << " not found!" << endl;
    exit(1);
  }

  // Store bin information

  map<string,size_t>columnNameMap;
  for (int i=0; i<*nBinDimension; i++) {
    columnNameMap[stringFromFortran(binNames[i],COLUMN_NAME_LEN)]=i;
  }
  TheorEval*te=it->second;
  te->setBins(*nBinDimension, *nPoints, binFlags, allBins,columnNameMap);
  return 1;
}

/*!
 Initializes theory for requested dataset.
 */
void init_theor_eval_(int *dsId)
{
  tTEmap::iterator it = gTEmap.find(*dsId);
  if (it == gTEmap.end() ) {
    cout << "ERROR: Theory evaluation for dataset ID " << *dsId
    << " not found!" << endl;
    exit(1);
  }

  TheorEval *te = gTEmap.at(*dsId);
  te->initTheory();
}

/*!
 Evaluates theory for requested dataset and writes it to the global THEO array.
 */
void get_theor_eval_(int *dsId, int *np, int*idx)
{

  tTEmap::iterator it = gTEmap.find(*dsId);
  if (it == gTEmap.end() ) {
    cout << "ERROR: Theory evaluation for dataset ID " << *dsId
    << " not found!" << endl;
    exit(1);
  }

  TheorEval *te = gTEmap.at(*dsId);
  valarray<double>vte(te->getNbins());//vector of theory predictions for this dataset
  te->Evaluate(vte);//writes into vte

  // write the predictions to THEO array
  const vector<int>*te_binflags=te->getBinFlags();
  const int*binflags=te_binflags->data();//get pointer to array of bin flags
  size_t ip=0;
  size_t offset=*idx-1;
  size_t endi=te_binflags->size();
  for(size_t i=0;i<endi;++i){
    if(binflags[i]!=0){//skip bins flagged 0
      c_theo_.theo[ip+offset]=vte[i];
      ++ip;
    }
  }

  if( ip != *np ){
    cout << "ERROR in get_theor_eval_: number of points mismatch" << endl;
    return;
  }
}

void close_theor_eval_()
{
  tTEmap::iterator it = gTEmap.begin();
  for (; it!= gTEmap.end(); it++){
    delete it->second;
  }

  gTEmap.clear();
}

void init_at_iteration_() {
  xfitter::updateDependentParameters();

  for(const auto it:XFITTER_PARS::gParameterisations){
    it.second->atIteration();
  }

  for ( auto pdfdecomposition : XFITTER_PARS::gPdfDecompositions) {
    pdfdecomposition.second->atIteration();//Among other things, sumrules are handled here
  }

  for(const auto it:XFITTER_PARS::gEvolutions){
    it.second->atIteration();
  }

  for(const auto it:XFITTER_PARS::gEvolutions){
    it.second->afterIteration();
  }

  for(const auto reaction:gNameReaction){
    reaction.second->atIteration();
  }
}
//This is called after minimization, after result output
//Could be named atEnd or something --Ivan
void fcn3action_()
{
  // Minimizer action:
  if (XFITTER_PARS::gMinimizer != nullptr ) {
    XFITTER_PARS::gMinimizer->actionAtFCN3();
  }

  for ( auto reaction : gNameReaction ) {
    reaction.second->atFCN3();
  }

  if (XFITTER_PARS::rootNode["ExtraActions"]["PrintPionMoments"].IsDefined()){
    {
    auto it = XFITTER_PARS::gParameters.find("Av");
    if (it != XFITTER_PARS::gParameters.end() ) {
      cout<<"  Av="<<*it->second<<endl;
    }
    it = XFITTER_PARS::gParameters.find("As");
    if (it != XFITTER_PARS::gParameters.end() ) {
      cout<<"  As="<<*it->second<<endl;
    }
    it = XFITTER_PARS::gParameters.find("Ag");
    if (it != XFITTER_PARS::gParameters.end() ) {
      cout<<"  Ag="<<*it->second<<endl;
    }
    }
    {
    auto it = XFITTER_PARS::gParameterisations.find("v");
    if (it != XFITTER_PARS::gParameterisations.end() ) {
      cout<<" <v>="<<it->second->moment(-1)<<endl;
      cout<<"<xv>="<<it->second->moment(0)<<endl;
    }
    it = XFITTER_PARS::gParameterisations.find("S");
    if (it != XFITTER_PARS::gParameterisations.end() ) {
      cout<<"<xS>="<<it->second->moment(0)<<endl;
    }
    it = XFITTER_PARS::gParameterisations.find("g");
    if (it != XFITTER_PARS::gParameterisations.end() ) {
      cout<<"<xg>="<<it->second->moment(0)<<endl;
    }
    }
  }
}

void error_band_action_(const int& i) {
  for ( auto reaction : gNameReaction ) {
    reaction.second->atIteration();
    reaction.second->atMakeErrorBands(i);
  }
}
