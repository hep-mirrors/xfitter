/*!
 @file theor_eval.cc
 @date Tue Aug 6 2013
 @author Andrey Sapronov <Andrey.Sapronov@cern.ch>

 Contains Fortran interface functions to operate with theoretical
 predictions obtained via fast cross section evaluation methods, 
 e.g. APPLgrid,  FastNLO and k-Factors.
 */

#include <vector>
#include <valarray>

#include "xfitter_cpp.h"

#include "TheorEval.h"
//#include "datasets.icc"

using namespace std;

extern "C" {
  int set_theor_eval_(int *dsId);//, int *nTerms, char **TermName, char **TermType, 
//    char **TermSource, char *TermExpr);
  int set_theor_bins_(int *dsId, int *nBinDimension, int *nPoints, int *binFlags, 
    double *allBins);
//  int set_theor_units_(int *dsId, double *units);
  int init_theor_eval_(int *dsId);
  int update_theor_ckm_();
  int get_theor_eval_(int *dsId, int* np, int* idx);
  int close_theor_eval_();
}

/// global dataset to theory evaluation pointer map
tTEmap gTEmap;


extern struct thexpr_cb {
  double dynscale;
  int nterms;
  char termname[16][8];
  char termtype[16][80];
  char terminfo[16][80];
  char termsource[16][1000];
  char theorexpr[1000];
  int ppbar_collisions;
  int normalised;
  int murdef;
  int mufdef;
} theorexpr_;

extern struct ord_scales {
   double datasetmur[150];
   double datasetmuf[150];
   int datasetiorder[150];
} cscales_;

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
  double *allBins)
{
  tTEmap::iterator it = gTEmap.find(*dsId);
  if (it == gTEmap.end() ) { 
    cout << "ERROR: Theory evaluation for dataset ID " << *dsId 
    << " not found!" << endl;
    exit(1);
  }
  
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
 Updates the CKM matrix to all the initialized appl grids
 */
int update_theor_ckm_()
{
  double a_ckm[] = { ckm_matrix_.Vud, ckm_matrix_.Vus, ckm_matrix_.Vub,
                                  ckm_matrix_.Vcd, ckm_matrix_.Vcs, ckm_matrix_.Vcb,
                                  ckm_matrix_.Vtd, ckm_matrix_.Vts, ckm_matrix_.Vtb};
  vector<double> v_ckm (a_ckm, a_ckm+sizeof(a_ckm)/sizeof(a_ckm[0]));
  tTEmap::iterator it = gTEmap.begin();
  for (; it!= gTEmap.end(); it++){
    it->second->setCKM(v_ckm);
  }
  
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
      c_theo_.theo_[*idx+ip-1]=vte[int(ibf-binflags->begin())];
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
