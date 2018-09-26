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

#include "get_pdfs.h"
#include "xfitter_cpp.h"

#include "TheorEval.h"
//#include "datasets.icc"
#include <yaml-cpp/yaml.h>
#include "ReactionTheory.h"
#include "xfitter_pars.h"

#include "BaseEvolution.h"
#include "BasePdfDecomposition.h"
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
  int init_theor_eval_(int *dsId);
  int update_theor_ckm_();
  int get_theor_eval_(int *dsId, int* np, int* idx);
  int read_reactions_();
  int close_theor_eval_();
  void init_func_map_();
  void init_at_iteration_(); ///< Loop over reactions, initialize them
  void fcn3action_();      ///< Loop over reactions, call actionAtFCN3
  void error_band_action_(const int& i); ///< Loop over rections, call error_band_action
}

/// global dataset to theory evaluation pointer map
tTEmap gTEmap;
tReactionLibsmap gReactionLibs;
tNameReactionmap gNameReaction;
tDataBins gDataBins;

t2Dfunctions g2Dfunctions;

extern struct thexpr_cb {
  double dynscale;
  int nterms;
  char termname[16][8];
  char termtype[16][80];
  char terminfo[16][2048];
  char termsource[16][256];
  char theorexpr[1000];
  int ppbar_collisions;
  int normalised;
  int murdef;
  int mufdef;
  int ninfo;  // dataset info as well  
  double datainfo[100];
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
int set_theor_eval_(int *dsId)//, int *nTerms, char **TermName, char **TermType, 
//  char **TermSource, char *TermExpr)
{
  // convert fortran strings to c++
  vector<string> stn(theorexpr_.nterms);
  vector<string> stt(theorexpr_.nterms);
  vector<string> sti(theorexpr_.nterms);
  vector<string> sts(theorexpr_.nterms);
  for ( int i = 0; i< theorexpr_.nterms; i++){
    stn[i].assign(theorexpr_.termname[i], string(theorexpr_.termname[i]).find(' '));
    stt[i].assign(theorexpr_.termtype[i], string(theorexpr_.termtype[i]).find(' '));
    sti[i].assign(theorexpr_.terminfo[i], string(theorexpr_.terminfo[i]).find(' '));
    sts[i].assign(theorexpr_.termsource[i], string(theorexpr_.termsource[i]).find(' '));
  }
  string ste;
  ste.assign(theorexpr_.theorexpr, string(theorexpr_.theorexpr).find(' '));
  TheorEval *te = new TheorEval(*dsId, theorexpr_.nterms, stn, stt, sti, sts, ste);

  // Store CINFO
  for (int i=0; i<theorexpr_.ninfo; i++) {
    std::string n(theorexpr_.CInfo[i]);
    n = n.substr(0,80);
    n.erase(std::remove(n.begin(), n.end(), ' '), n.end());
    te->AddDSParameter(n, theorexpr_.datainfo[i]);
  } 
  // Store some other basic info
  theorexpr_.dsname[79] = '\0';
  std::string n(theorexpr_.dsname);
  // Erase trailing spaces
  n = rtrim(n)," ";

  te->SetDSname(n);
  te->AddDSParameter("Index",theorexpr_.ds_index); // dataset index
  te->AddDSParameter("FileIndex",*dsId); 

  te->SetCollisions(theorexpr_.ppbar_collisions);
  te->SetDynamicScale(theorexpr_.dynscale);
  te->SetNormalised(theorexpr_.normalised);
  te->SetMurMufDef(theorexpr_.murdef,theorexpr_.mufdef);
  te->SetOrdScales(cscales_.datasetiorder[*dsId-1],cscales_.datasetmur[*dsId-1],cscales_.datasetmuf[*dsId-1]);

  tTEmap::iterator it = gTEmap.find(*dsId);
  if (it == gTEmap.end() ) { gTEmap[*dsId] = te; }
  else {
    cout << "ERROR: Theory evaluation for dataset ID " << *dsId 
    << " already exists." << endl;
    exit(1); // make proper exit later
  }

  return 1;
}

/*!
 Sets datasets bins in theory evaluations.
 write details on argumets
 */
int set_theor_bins_(int *dsId, int *nBinDimension, int *nPoints, int *binFlags, 
		    double *allBins, char binNames[10][80])
{
  tTEmap::iterator it = gTEmap.find(*dsId);
  if (it == gTEmap.end() ) { 
    cout << "ERROR: Theory evaluation for dataset ID " << *dsId 
    << " not found!" << endl;
    exit(1);
  }
  
  // Store bin information

  map<string, valarray<double> > namedBins;
  for (int i=0; i<*nBinDimension; i++) {
    string name = binNames[i];
    name.erase(name.find(" "));
    //    cout << name << " " << *dsId <<endl;
    valarray<double> bins(*nPoints); 
    for ( int j = 0; j<*nPoints; j++) {
      bins[j] = allBins[j*10 + i];
    }

    //namedBins[name] = bins; // OZ 30.012017 this is not legal in C++ < C++11 and does not work with gcc 4.4.7
    valarray<double>& insertedBins = namedBins[name];
    insertedBins.resize(bins.size());
    insertedBins = bins;
  }
  gDataBins[*dsId] = namedBins;

  TheorEval *te = gTEmap.at(*dsId);
  te->setBins(*nBinDimension, *nPoints, binFlags, allBins);
  return 1;
}

/*
int set_theor_units_(int *dsId, double *units)
{
  tTEmap::iterator it = gTEmap.find(*dsId);
  if (it == gTEmap.end() ) { 
    cout << "ERROR: Theory evaluation for dataset ID " << *dsId 
    << " not found!" << endl;
    exit(1);
  }
  
  TheorEval *te = gTEmap.at(*dsId);
  te->setUnits(*units);
  return 1;
}
*/

/*!
 Initializes theory for requested dataset.
 */
int init_theor_eval_(int *dsId)
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
int get_theor_eval_(int *dsId, int *np, int*idx)
{

  tTEmap::iterator it = gTEmap.find(*dsId);
  if (it == gTEmap.end() ) { 
    cout << "ERROR: Theory evaluation for dataset ID " << *dsId 
    << " not found!" << endl;
    exit(1);
  }
  
  valarray<double> vte;
  TheorEval *te = gTEmap.at(*dsId);
  vte.resize(te->getNbins());
  te->Evaluate(vte);

  // Get bin flags, and abandon bins flagged 0
  const vector<int> *binflags = te->getBinFlags();
  int ip = 0;
  vector<int>::const_iterator ibf = binflags->begin();
  for (; ibf!=binflags->end(); ibf++){
    if ( 0 != *ibf ) {
      c_theo_.theo[*idx+ip-1]=vte[int(ibf-binflags->begin())];
      ip++;
    }
      //cout << *ibf << "\t" << vte[int(ibf-binflags->begin())] << endl;
  }

  // write the predictions to THEO array
  if( ip != *np ){
    cout << "ERROR in get_theor_eval_: number of points mismatch" << endl;
    return -1;
  }
}

int close_theor_eval_()
{
  tTEmap::iterator it = gTEmap.begin();
  for (; it!= gTEmap.end(); it++){
    delete it->second;
  }

  gTEmap.clear();
}


/*!
 */
int read_reactions_()
{
  ifstream frt((PREFIX+string("/lib/Reactions.txt")).c_str());
  if ( frt.is_open() ) {
    while (1){
      string rname, lib;
      frt >> rname >> lib;
      if (frt.eof()) break;
      if (gReactionLibs.find(rname) == gReactionLibs.end() ) {
	// possible check
      }
      gReactionLibs[rname] = lib;

    }
  }
  else {
    string text = string("F: can not open Reactions.txt file. Check your ")+PREFIX+string("/lib directory");
    hf_errlog_(16121401,text.c_str(),text.size());
  }
  return 1;
}


// a bunch of functions 
double xg(const double& x, const double& q2) {  double pdfs[20]; HF_GET_PDFS_WRAP(x,q2,pdfs); return pdfs[6+0]; }
double xu(const double& x, const double& q2) {  double pdfs[20]; HF_GET_PDFS_WRAP(x,q2,pdfs); return pdfs[6+1]; }
double xub(const double& x, const double& q2) {  double pdfs[20]; HF_GET_PDFS_WRAP(x,q2,pdfs); return pdfs[6-1]; }


void init_func_map_() {
  g2Dfunctions["xg"] = &xg;
  g2Dfunctions["xu"] = &xu;
  g2Dfunctions["xub"] = &xub;
}

void init_at_iteration_() {
  
  for ( auto pdfdecomposition : XFITTER_PARS::gPdfDecompositions) {
    pdfdecomposition.second->atIteration();
  }

  
  for(auto it:XFITTER_PARS::gEvolutions) {
    xfitter::BaseEvolution*evolution=it.second;
    evolution->atIteration();

    // register updated PDF XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    //Wait, do they even get updated between iterations? Is this here even necessary? --Ivan
    
    XFITTER_PARS::registerXfxQArray(evolution->_name,evolution->xfxQArray());
  }


 for ( auto reaction : gNameReaction ) {
    reaction.second->initAtIteration();
  }

}

void fcn3action_()
{
  // Minimizer action:
  if (XFITTER_PARS::gMinimizer != nullptr ) {
    XFITTER_PARS::gMinimizer->actionAtFCN3();
  }
  
  for ( auto reaction : gNameReaction ) {
    reaction.second->actionAtFCN3();
  }
}

void error_band_action_(const int& i) {
  for ( auto reaction : gNameReaction ) {
    reaction.second->initAtIteration();   // Init parameters first
    reaction.second->errorBandAction(i);
  }
}
