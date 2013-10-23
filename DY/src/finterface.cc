#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

#include "PhysPar.h"
#include "DYcalc.h"
#include "IntSteps.h"
#include "PDFconv.h"
#include "BinMatrix.h"
#include <stdlib.h>

using namespace std;

extern "C" {
  int dy_create_calc_(const int *ds_id, const int *chg_prod, 
      const double *beam_en, const char *boz,
      const double *ranges, const char *var_name, 
      const int *n_bins, const double *bin_edges);

  int dy_do_calc_();

  int dy_get_res_(const int *ds_id, double *calc_res);

  int dy_release_();
  int dy_set_ewpars_();
}

typedef map <int, DYcalc* > DCmap;
DCmap     gCalcs;
vector<PDFconv*> gPDFconvs;
vector<BinMatrix*> gBinMatrices;

// initializes Drell-Yan LO calculations with info on
// beam, process, kinematic cuts, and bins.
int dy_create_calc_(const int *ds_id, const int *chg_prod, 
    const double *beam_en, const char *boz,
    const double *ranges, const char *var_name, 
    const int *n_bins, const double *bin_edges)
{

  // initialize integration grid if didn't yet.
  const IntSteps *int_steps = new IntSteps(string(boz), ranges, 
    var_name, *n_bins, bin_edges);

  BinMatrix *bm = NULL;
  // check for existing bin matrices
  DCmap::iterator idc = gCalcs.begin(); 

//  for ( ;idc!=gCalcs.end();idc++ ) {
//    BinMatrix *tbm = idc->second->getBM();
//    if ( (tbm->getBeamEn() == *beam_en )
//      && (string(boz) == tbm->getBozName()) ){
//      bm = tbm;
//      /*
//      if ( string(boz) != bm->getBozName() ) {
//        cout << "Simultaneous calculation of DY charged and neutral \n\
//	         is not supported. Abort. " << endl;
//	exit(1);
//      }
//      */
//      break;
//    }
//  }

  // if not found create new
  if ( NULL == bm ) {
    bm = new BinMatrix(beam_en, int_steps);
  }

  // prepare PDF grids
  PDFconv *pc = NULL;
  // check for existing bin matrices
  idc = gCalcs.begin(); 
//  for ( ;idc!=gCalcs.end();idc++ ) {
//    PDFconv *tpc = idc->second->getPC();
//    if ( tpc->isSameBeam(*chg_prod, beam_en) && 
//         (string(boz) == tpc->getBozName() )){
//      pc = tpc;
//      /*
//      if ( string(boz) != pc->getBozName() ) {
//        cout << "Simultaneous calculation of DY charged and neutral \n\
//	         is not supported. Abort. " << endl;
//	exit(1);
//      }
//      */
//      break;
//    }
//  }

  // if not found create new
  if ( NULL == pc ) {
    pc = new PDFconv(*chg_prod, beam_en, int_steps);
    gPDFconvs.push_back(pc);
  }


  // create calculator and put to map
  DYcalc * dc = new DYcalc(bm, pc, int_steps);
  gCalcs.insert( pair<int,DYcalc*>( *ds_id,dc ) );

  return 1;
}


// calculate Drell-Yan LO cross sections for all data sets
int dy_do_calc_()
{
  // evolve convolutions
  vector<PDFconv*>::iterator ipc = gPDFconvs.begin();
  for (; ipc!=gPDFconvs.end(); ipc++){
    (*ipc)->interpPDF();
  }

  DCmap::iterator idc = gCalcs.begin();
  for (; idc != gCalcs.end(); idc++ ){
    if ( true != idc->second->Integrate() ) {
      cout << "Something is wrong with DY integration for " 
           << idc->first << " data set." << endl;
      return 0;
    }
  }

  return 1;
}


// return DY calculations for data set ds_name
int dy_get_res_(const int *ds_id, double *calc_res)
{
  DYcalc * dc = gCalcs.find(*ds_id)->second;
  dc->getCalcRes(calc_res);

  return 1;
}

int dy_set_ewpars_(){
  PhysPar::setPhysPar();
}

int dy_release_()
{
  vector<PDFconv*>::iterator ipc = gPDFconvs.begin();
  for (; ipc!=gPDFconvs.end(); ipc++){
    delete *ipc;
  }
  
  vector<BinMatrix*>::iterator ibm = gBinMatrices.begin();
  for (; ibm!=gBinMatrices.end(); ibm++){
    delete *ibm;
  }

  DCmap::iterator idc = gCalcs.begin();
  for (; idc != gCalcs.end() ; idc++){
    delete (idc->second);
  }

  return 1;
}
