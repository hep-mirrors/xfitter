
/*
   @file ReactionHVQMNR_LHCb_7TeV_charm.cc
   @date 2017-01-02
   @author Oleksandr Zenaiev (oleksandr.zenaiev@desy.de)
   Created by  AddReaction.py on 2017-01-02

   Derived from ReactionBaseHVQMNR where basic stuff for HVQMNR calculation is implemented
   This class implements calculation for LHCb charm measurement at 7 TeV
   [Nucl. Phys. B 871 (2013), 1] [arXiv:1302.2864]
*/

#include "ReactionHVQMNR_LHCb_7TeV_charm.h"
#include "xfitter_cpp.h"
#include <TF1.h>
#include <TH2D.h>
#include <TMath.h>

// the class factories
extern "C" ReactionHVQMNR_LHCb_7TeV_charm* create() {
  return new ReactionHVQMNR_LHCb_7TeV_charm();
}


// initialize at the start of the computation
int ReactionHVQMNR_LHCb_7TeV_charm::atStart(const string &s)
{
  // ignore provided terminfo (s): all needed information has been set already
  // via setDatasetParameters(int dataSetID, map<string,string> pars)

  // ******************************************************************
  // perform initialisation and pre-calculation
  // ******************************************************************
  // protection against overdoing
  if(_isInitAtStart)
    return 0;
  _isInitAtStart = true;
  //printf("ReactionHVQMNR_LHCb_7TeV_charm::atStart()\n");

  // check HF scheme
  CheckHFScheme();

  // read needed theory parameters
  UpdateParameters();
  PrintParameters();

  // stereing parameters for this calculation (modify only if you understand what you are doing)
  Steering steer;
  steer.ptmin = 0.001;
  steer.ptmax = 20.0;
  steer.npt   = 25;
  steer.nptsm = 500;
  steer.ymin  = 2.0;
  steer.ymax  = 6.0;
  steer.ny    = 40;
  steer.nsfnb = 500;
  steer.nx3   = 25;
  steer.nx4   = 125;
  steer.nbz   = 50;

  DefaultInit(steer, _pars.mc, _mnr, _frag, _grid, _gridSmoothed);
  //if(_debug)
  //printf("ReactionHVQMNR_LHCb_7TeV_charm::atStart(): at initialisation mc = %f\n", _pars.mc);
  // MNR (parton-level calculation)
  _mnr.SetDebug(_debug);
  _mnr.fC_sh = TMath::Power(7000.0, 2.0); // centre-of-mass energy squared
  _mnr.CalcConstants();
  // fragmentation
  _frag.SetDebug(_debug);
  // add needed final states (data set specific)
  _frag.AddOut(MNR::Frag::GetFragFunction(0, "dzero", _pars.fragpar_c), MNR::Frag::fM_dzero);
  _frag.AddOut(MNR::Frag::GetFragFunction(0, "",      _pars.fragpar_c), MNR::Frag::fM_dstar);
  _frag.AddOut(MNR::Frag::GetFragFunction(0, "dch",   _pars.fragpar_c), MNR::Frag::fM_dch);
  _frag.AddOut(MNR::Frag::GetFragFunction(0, "",      _pars.fragpar_c), MNR::Frag::fM_ds);
  _frag.AddOut(MNR::Frag::GetFragFunction(0, "",      _pars.fragpar_c), MNR::Frag::fM_lambdac);
  _frag.AddOut(MNR::Frag::GetFragFunction(0, "",      _pars.fragpar_c), MNR::Frag::fM_lambdac);

  // Set binning for cross-section histograms (data set specific)
  _hCalculatedXSec.resize(6);
  for(int f = 0; f < 6; f++)
    _hCalculatedXSec[f] = new TH2D;
  // D0, D+, D*+, D_s
  int nbin_y = 5;
  double bin_y[6] = { 2.0, 2.5, 3.0, 3.5, 4.0, 4.5 };
  int nbin_pt = 8;
  double bin_pt[9] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
  for(int f = 0; f < 4; f++)
    _hCalculatedXSec[f]->SetBins(nbin_pt, bin_pt, nbin_y, bin_y);
  // Lambda_c rapidity differenential (for normalised cross section)
  int nbin_pt_lambdac = 1;
  double bin_pt_lambdac[2] = { 2.0, 8.0 };
  _hCalculatedXSec[4]->SetBins(nbin_pt_lambdac, bin_pt_lambdac, nbin_y, bin_y);
  // Lambda_c pT differenential (for unnormalised cross section)
  int nbin_y_lambdac = 1;
  double bin_y_lambdac[2] = { 2.0, 4.5 };
  _hCalculatedXSec[5]->SetBins(nbin_pt, bin_pt, nbin_y_lambdac, bin_y_lambdac);

  return 0;
}


// perform calculation (this is done once per iteration)
void ReactionHVQMNR_LHCb_7TeV_charm::initAtIteration()
{
  // protection against overdoing
  // TODO: remove this trick
  //if(_ifcncount_last == cfcn_.ifcncount)
  //  return;
  //_ifcncount_last = cfcn_.ifcncount;

  //printf("ReactionHVQMNR_LHCb_7TeV_charm::initAtIteration() %d\n", cfcn_.ifcncount);

  // read needed MINUIT parameters (enough to be done once per iteration)
  UpdateParameters();
  if(_debug)
    PrintParameters();

  // update parameters and perform calculation
  _mnr.SetScaleCoef(_pars.mf_A_c, _pars.mf_B_c, _pars.mf_C_c, _pars.mr_A_c, _pars.mr_B_c, _pars.mr_C_c);
  _mnr.CalcXS(&_grid, _pars.mc);
  MNR::Grid::InterpolateGrid(&_grid, &_gridSmoothed, _pars.mc);
  for(int f = 0; f < 6; f++)
    _frag.GetFF(f)->SetParameter(1, _pars.fragpar_c);
  _frag.CalcCS(&_gridSmoothed, _pars.mc, _hCalculatedXSec);
}


// main function to compute results at an iteration
int ReactionHVQMNR_LHCb_7TeV_charm::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
  // TODO move to core xFitter
  //initAtIteration();
  //printf("ReactionHVQMNR_LHCb_7TeV_charm::compute() %d\n", dataSetID);

  // get histogramm with cross sections for needed dataset
  DataSet& ds = _dataSets[dataSetID];
  TH2D* histXSec = NULL;
  if(ds.FinalState == "dzero")
    histXSec = _hCalculatedXSec[0];
  else if(ds.FinalState == "dstar")
    histXSec = _hCalculatedXSec[1];
  else if(ds.FinalState == "dch")
    histXSec = _hCalculatedXSec[2];
  else if(ds.FinalState == "ds")
    histXSec = _hCalculatedXSec[3];
  else if(ds.FinalState == "lambdac" && ds.NormY == 1)
    histXSec = _hCalculatedXSec[4];
  else if(ds.FinalState == "lambdac" && ds.NormY == 0)
    histXSec = _hCalculatedXSec[5];
  else
    hf_errlog(16123005, "F: in HVQMNR_LHCb_7tev_charm(): unknown FinalState" + ds.FinalState);

  // match bins and fill cross section array
  bool diffPt = false;
  bool diffY = false;
  for(unsigned int i = 0; i < ds.BinsYMin->size(); i++)
  {
    val[i] = FindXSecPtYBin(histXSec, (*ds.BinsYMin)[i], (*ds.BinsYMax)[i], (*ds.BinsPtMin)[i], (*ds.BinsPtMax)[i], diffPt, diffY);
    // normalise to another bin if needed
    if(ds.NormY == 1)
      val[i] = val[i] / FindXSecPtYBin(histXSec, (*ds.BinsYMinRef)[i], (*ds.BinsYMaxRef)[i], (*ds.BinsPtMin)[i], (*ds.BinsPtMax)[i], diffPt, diffY);
    // ... or scale to fragmentation fraction
    else
      val[i] = val[i] * ds.FragFraction;
  }

  return 0;
}
