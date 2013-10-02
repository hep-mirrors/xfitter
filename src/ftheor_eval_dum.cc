/*!
 @file theor_eval.cc
 @date Tue Aug 6 2013
 @author Andrey Sapronov <Andrey.Sapronov@cern.ch>

 Contains Fortran interface functions to operate with theoretical
 predictions obtained via fast cross section evaluation methods, 
 e.g. APPLgrid,  FastNLO and k-Factors.
 */

#include <iostream>

using namespace std;

extern "C" {
  int set_theor_eval_(int *dsId);//, int *nTerms, char **TermName, char **TermType, 
//    char **TermSource, char *TermExpr);
  int set_theor_bins_(int *dsId, int *nBinDimension, int *nPoints, int *binFlags, 
    double *allBins);
  int set_theor_units_(int *dsId, double *units);
  int init_theor_eval_(int *dsId);
  int get_theor_eval_(int *dsId, int *iorder, double *mur, double *muf, 
    int* np, int* idx);
  int close_theor_eval_();
}

/*!
 Creates theory evaluation object and adds it to the global map by 
 dataset ID.
 write details on argumets
 */
int set_theor_eval_(int *dsId)//, int *nTerms, char **TermName, char **TermType, 
//  char **TermSource, char *TermExpr)
{
  cout << "ERROR: calling dummy procedure for theory evaluation in dataset " << *dsId 
  << ". Recompile with --enable-applgrid to use this option." << endl;

  return -1;
}

/*!
 Sets datasets bins in theory evaluations.
 write details on argumets
 */
int set_theor_bins_(int *dsId, int *nBinDimension, int *nPoints, int *binFlags, 
  double *allBins)
{
  cout << "ERROR: calling dummy procedure for theory evaluation in dataset " << *dsId 
  << ". Recompile with --enable-applgrid to use this option." << endl;

  return -1;
}

int set_theor_units_(int *dsId, double *units)
{
  cout << "ERROR: calling dummy procedure for theory evaluation in dataset " << *dsId 
  << ". Recompile with --enable-applgrid to use this option." << endl;

  return -1;
}

/*!
 Initializes theory for requested dataset.
 */
int init_theor_eval_(int *dsId)
{
  cout << "ERROR: calling dummy procedure for theory evaluation in dataset " << *dsId 
  << ". Recompile with --enable-applgrid to use this option." << endl;

  return -1;
}

/*!
 Evaluates theory for requested dataset and writes it to the global THEO array.
 */
int get_theor_eval_(int *dsId, int *iorder, double *mur, double *muf, int *np, int*idx)
{
  cout << "ERROR: calling dummy procedure for theory evaluation in dataset " << *dsId 
  << ". Recompile with --enable-applgrid to use this option." << endl;

  return -1;
}

int close_theor_eval_()
{
  return 1;
}
