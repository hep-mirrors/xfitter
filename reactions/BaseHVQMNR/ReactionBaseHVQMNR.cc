#include "xfitter_cpp.h"
#include <TMath.h>
#include <string>

/*
   @file ReactionBaseHVQMNR.cc
   @date 2017-01-02
   @author Oleksandr Zenaiev (oleksandr.zenaiev@desy.de)
   Created by  AddReaction.py on 2017-01-02
   *
   This is abstract class from which implementations of HVQMNR
   calculations for particular datasets should be derived
*/

#include "ReactionBaseHVQMNR.h"

// the class factories
extern "C" ReactionBaseHVQMNR* create() {
  // this is abstract class, no instance should be created
  //return new ReactionBaseHVQMNR();
  hf_errlog(16123000, "F: can not create ReactionBaseHVQMNR instance: you should implement calculation in a derived class");
  return NULL;
}


// pass to MNR pointer to instance inherited from ReactionTheory to allow access to alphas and PDF routines
ReactionBaseHVQMNR::ReactionBaseHVQMNR() : _mnr(MNR::MNR(this))
{
  // set initialisation status flag
  _isInitAtStart = false;
  // set debugging flag
  _debug = steering_.ldebug;
  //_debug = 1;
}


ReactionBaseHVQMNR::~ReactionBaseHVQMNR()
{
  for(unsigned int i = 0; i < _hCalculatedXSec.size(); i++)
    delete _hCalculatedXSec[i];
}


void ReactionBaseHVQMNR::initTerm(TermData *td)
{
  unsigned dataSetID = td->id;
  _tdDS[dataSetID] = td;

  // Order
  //const string order = td->getParamS("Order");
  //if(order != "NLO")
  //  hf_errlog(19020301, "F: order " + order + " not supported");

  // add new dataset
  DataSet dataSet;
  std::pair<std::map<int, DataSet>::iterator, bool> ret = _dataSets.insert(std::pair<int, DataSet>(dataSetID, dataSet));
  // check if dataset with provided ID already exists
  if(!ret.second)
    hf_errlog(16123001, "F: dataset with id = " + std::to_string(dataSetID) + " already exists");

  // set parameters for new dataset
  DataSet& ds = ret.first->second;
  // mandatory "FinalState="
  if(!td->hasParam("FinalState"))
    hf_errlog(16123002, "F: TermInfo must contain FinalState entry");
  else
    ds.FinalState = td->getParamS("FinalState");

  // optional "NormY="
  if(!td->hasParam("NormY"))
    ds.NormY = 0; // default value is unnormalised absolute cross section
  else
    ds.NormY = td->getParamI("NormY");

  // optional "FragFrac=" (must be provided for absolute cross section)
  if(!td->hasParam("FragFrac"))
  {
    if(ds.NormY == 0)
      hf_errlog(16123003, "F: for absolute cross section TermInfo must contain FragFrac entry");
  }
  else
  {
    ds.FragFraction = *td->getParamD("FragFrac");
    if(ds.NormY != 0)
      printf("Warning: FragFrac=%f will be ignored for normalised cross sections\n", ds.FragFraction);
  }

  // set binning
  ds.BinsYMin  = const_cast<std::valarray<double>*>(td->getBinColumnOrNull("ymin"));
  ds.BinsYMax  = const_cast<std::valarray<double>*>(td->getBinColumnOrNull("ymax"));
  ds.BinsPtMin = const_cast<std::valarray<double>*>(td->getBinColumnOrNull("pTmin"));
  ds.BinsPtMax = const_cast<std::valarray<double>*>(td->getBinColumnOrNull("pTmax"));
  //if (ds.BinsYMin == NULL || ds.BinsYMax == NULL || ds.BinsPtMin == NULL || ds.BinsPtMax == NULL )
  //  hf_errlog(16123004, "F: No bins ymin or ymax or ptmin or ptmax");
  // set reference y bins if needed
  if(ds.NormY == 1)
  {
    ds.BinsYMinRef = const_cast<std::valarray<double>*>(td->getBinColumnOrNull("yminREF"));
    ds.BinsYMaxRef = const_cast<std::valarray<double>*>(td->getBinColumnOrNull("ymaxREF"));
    if(ds.BinsYMinRef == NULL || ds.BinsYMaxRef == NULL)
      hf_errlog(16123005, "F: No bins yminREF or ymaxREF for normalised cross section");
  }

  if(_debug)
    printf("Added dataset: FinalState = %s  NormY = %d  FragFraction = %f\n", ds.FinalState.c_str(), ds.NormY, ds.FragFraction);
}


// ********************************
// ***** utility routines *********
// ********************************

// check equality of float numbers with tolerance
bool ReactionBaseHVQMNR::IsEqual(const double val1, const double val2, const double eps/* = 1e-6*/)
{
  return (TMath::Abs(val1 - val2) < eps || TMath::Abs((val1 - val2) / val1) < eps);
}

// initialise calculation with default parameters
void ReactionBaseHVQMNR::DefaultInit(const Steering& steer, const double mq, MNR::MNR& mnr, MNR::Frag& frag, MNR::Grid& grid, MNR::Grid& grid_smoothed)
{
  DefaultInitMNR(steer, mq, mnr);
  DefaultInitGrid(steer, mq, steer.npt, grid);
  DefaultInitGrid(steer, mq, steer.nptsm, grid_smoothed);
  DefaultInitFrag(steer, frag);
}

void ReactionBaseHVQMNR::DefaultInitMNR(const ReactionBaseHVQMNR::Steering &steer, const double mq, MNR::MNR &mnr)
{
  // MNR parton level cross sections, quark-antiquark contributions
  mnr.bFS_Q = steer.q;
  mnr.bFS_A = steer.a;
  // number of light flavours
  mnr.fC_nl = steer.nf;
  // x3 and x4 binning
  mnr.fBn_x3 = steer.nx3;
  mnr.fBn_x4 = steer.nx4;
  mnr.fSF_nb = steer.nsfnb;
  // PDF range
  mnr.fSF_min_x = steer.xmin;
  mnr.fSF_max_x = steer.xmax;
  mnr.fSF_min_mf2 = steer.mf2min;
  mnr.fSF_max_mf2 = steer.mf2max;
  // precalculation (memory allocation etc.)
  mnr.CalcBinning();
}

void ReactionBaseHVQMNR::DefaultInitGrid(const ReactionBaseHVQMNR::Steering &steer, const double mq, const int npt, MNR::Grid &grid)
{
  // Parton level pT-y grids
  grid.SetL(npt, steer.ptmin, steer.ptmax, mq);
  grid.SetY(steer.ny, steer.ymin, steer.ymax);
  grid.SetW(1);
}

void ReactionBaseHVQMNR::DefaultInitFrag(const ReactionBaseHVQMNR::Steering &steer, MNR::Frag &frag)
{
  // Fragmentation
  frag.SetNz(steer.nbz);
}

// return cross section in provided pT-y bin
double ReactionBaseHVQMNR::FindXSecPtYBin(const TH2* histXSec, const double ymin, const double ymax, const double ptmin, const double ptmax, const bool diff_pt, const bool diff_y)
{
  // **************************************************************************
  // Input parameters:
  //   histXSec: 2D historam with calculated cross section (x axis pT, y axis rapidity)
  //   ymin: min rapidity of required bin
  //   ymax: max rapidity of required bin
  //   ptmin: min pT of required bin
  //   ptmax: max pT of required bin
  //   d_pt: if TRUE, cross section is divided by pT bin width
  //   d_y: if TRUE, cross section is divided by rapidity bin width
  // **************************************************************************
  //printf("FindXSecPtYBin: y %f %f pT %f %f\n", ymin, ymax, ptmin, ptmax);
  //histXSec->Print("all");
  for(int bpt = 0; bpt < histXSec->GetNbinsX(); bpt++)
  {
    //printf("bpt = %d\n", bpt);
    for(int by = 0; by < histXSec->GetNbinsY(); by++)
    {
      //printf("by = %d\n", by);
      if(!IsEqual(ymin, histXSec->GetYaxis()->GetBinLowEdge(by + 1))) continue;
      if(!IsEqual(ymax, histXSec->GetYaxis()->GetBinUpEdge(by + 1))) continue;
      if(!IsEqual(ptmin, histXSec->GetXaxis()->GetBinLowEdge(bpt + 1))) continue;
      if(!IsEqual(ptmax, histXSec->GetXaxis()->GetBinUpEdge(bpt + 1))) continue;
      double val = histXSec->GetBinContent(bpt + 1, by + 1);
      // check if XS is not nan
      if(val != val)
      {
        printf("Warning: nan cross section in bin y: %f %f, pt: %f %f, resetting to -1000000000.0\n", ymin, ymax, ptmin, ptmax);
        val = -1000000000.0;
      }
      // divide by bin width if needed
      if(diff_pt)
        val = val / (ptmax - ptmin);
      if(diff_y)
        val = val / (ymax - ymin);
      return val;
    }
  }
  hf_errlog(16123006, "F: ERROR in FindXSPtY(): bin not found y: " + std::to_string(ymin) + " " + std::to_string(ymax) +
            " , pt: " + std::to_string(ptmin) + " " + std::to_string(ptmax));
  return -1.0;
}

// check if appropriate heavy-flavour scheme is used
void ReactionBaseHVQMNR::CheckHFScheme()
{
  // check HF scheme
  //if(steering_.hfscheme != 3 && steering_.hfscheme != 4)
  //  hf_errlog(16123007, "S: calculation does not support HFSCHEME = " + std::to_string(steering_.hfscheme) + " (only 3, 4 supported)");
}

// read parameters for perturbative scales from MINUIT extra parameters
void ReactionBaseHVQMNR::GetMuPar(TermData* td, const char mu, const char q, double& A, double& B, double& C)
{
  // ***********************************************************************************************
  // Scales for charm and beauty production are parametrised as:
  // mu_f(c)^2 = MNRmf_A_c * pT_c^2 + MNRmf_B_c * m_c^2 + MNRmf_C_c
  // mu_r(c)^2 = MNRmr_A_c * pT_c^2 + MNRmr_B_c * m_c^2 + MNRmr_C_c
  // mu_f(b)^2 = MNRmf_A_b * pT_b^2 + MNRmf_B_b * m_b^2 + MNRmf_C_b
  // mu_r(b)^2 = MNRmr_A_b * pT_b^2 + MNRmr_B_b * m_b^2 + MNRmr_C_b
  // where mu_f(c), mu_r(c), mu_f(b), mu_r(b) are factorisation and renormalisation
  // scales for charm and beauty production, respectively, pT is transverse momentum
  // and m_c, m_b are charm and beauty quark masses.
  //
  // In total, one can provide all 12 parameters (MNRmf_A_c, MNRmf_B_c, MNRmf_C_c,
  // MNRmr_A_c, MNRmr_B_c, MNRmr_C_c, MNRmf_A_b, MNRmf_B_b, MNRmf_C_b,
  // MNRmr_A_b, MNRmr_B_b, MNRmr_C_b), however there are the foolowing rules:
  // 1) if suffix _c (_b) at the end of variable name is omitted, the parameter is applied
  //    for both charm and beauty production
  // 2) instead of providing e.g. MNRmr_A_c and MNRmr_B_c the user can provide one
  //    parameter MNRmr_AB_c so then MNRmr_A_c = MNRmr_B_c = MNRmr_AB_c
  // 3) if parameters *_C_* are not provided, they are set to 0
  // So e.g. for charm production only it is enough to provide just the following two parameters
  // MNRmr_AB and MNRmf_AB.
  // input variables:
  //   mu should be either 'f' or 'r' (for factorisation or renormalisation scale, respectively)
  //   q should be either 'c' or 'b' (for charm or beauty, respectively)
  // ***********************************************************************************************

  // ***************************
  const double defA = 1.0;
  const double defB = 1.0;
  const double defC = 0.0;
  // ***************************
  std::string baseParameterName = "MNRm" + std::string(1, mu);

  // A and B parameters
  if(td->hasParam(baseParameterName + "_AB"))
    A = B = *td->getParamD(baseParameterName + "_AB");
  else
  {
    if(td->hasParam(baseParameterName + "_A") && td->hasParam(baseParameterName + "_B"))
    {
      A = *td->getParamD(baseParameterName + "_A");
      B = *td->getParamD(baseParameterName + "_B");
    }
    else
    {
      if(td->hasParam(baseParameterName + "_AB_" + std::string(1, q)))
        A = B = *td->getParamD(baseParameterName + "_AB_" + std::string(1, q));
      else
      {
        if(td->hasParam(baseParameterName + "_A_" + std::string(1, q)))
          A = *td->getParamD(baseParameterName + "_A_" + std::string(1, q));
        else
          A = defA;
        if(td->hasParam(baseParameterName + "_B_" + std::string(1, q)))
          B = *td->getParamD(baseParameterName + "_B_" + std::string(1, q));
        else
          B = defB;
      }
    }
  }
  // C parameter
  if(td->hasParam(baseParameterName + "_C"))
    C = *td->getParamD(baseParameterName + "_C");
  else
  {
    if(td->hasParam(baseParameterName + "_C_" + std::string(1, q)))
      C = *td->getParamD(baseParameterName + "_C_" + std::string(1, q));
    else
      C = defC;
  }
}


// read fragmentation parameter from MINUIT extra parameters
double ReactionBaseHVQMNR::GetFragPar(TermData* td, const char q, const map<string,string> pars)
{
  // *********************************************************************
  // Parameters for non-perturbative fragmentation can be provided
  // as MINUIT extra parameters MNRfrag_c or MNRfrag_b
  // for charm and beauty production, respectively.
  // q should be either 'c' or 'b' (for charm or beauty, respectively)
  // *********************************************************************

  // ***************************
  const double defFFc = 4.4;
  const double defFFb = 11.0;
  // ***************************
  double parvalue = NAN;
  char parname[16];
  sprintf(parname, "MNRfrag_%c", q);
  if(!td->hasParam(parname))
  {
    // parameter not in ExtraParamMinuit -> using default value
    if(q == 'c')
      parvalue = defFFc;
    else if(q == 'b')
      parvalue = defFFb;
    else
      hf_errlog(17102103, "F: no default value for q = " + std::string(1, q) + " in ReactionBaseHVQMNR::GetFragPar()");
  }
  else
    parvalue = *td->getParamD(parname);
  return parvalue;
}

// read and update theory parameters
void ReactionBaseHVQMNR::UpdateParameters(TermData *td)
{
  // if not TermData provided, take any TermData pointer to access theory parameters
  // (theory parameters are supposed to be universal for all data sets in this case)
  if(!td)
    td = _tdDS.begin()->second;

  // heavy-quark masses
  _pars.mc = *td->getParamD("mch");
  _pars.mb = *td->getParamD("mbt");
  // scale parameters
  GetMuPar(td, 'f', 'c', _pars.mf_A_c, _pars.mf_B_c, _pars.mf_C_c);
  GetMuPar(td, 'r', 'c', _pars.mr_A_c, _pars.mr_B_c, _pars.mr_C_c);
  GetMuPar(td, 'f', 'b', _pars.mf_A_b, _pars.mf_B_b, _pars.mf_C_b);
  GetMuPar(td, 'r', 'b', _pars.mr_A_b, _pars.mr_B_b, _pars.mr_C_b);
  // fragmentation parameters
  _pars.fragpar_c = GetFragPar(td, 'c');
  _pars.fragpar_b = GetFragPar(td, 'b');

  // protection against not positive or nan masses
  if(_pars.mc <= 0.0 || _pars.mc != _pars.mc)
    _pars.mc = 1000.0;
  if(_pars.mb <= 0.0 || _pars.mb != _pars.mb)
    _pars.mb = 1000.0;
}

// print theory parameters
void ReactionBaseHVQMNR::PrintParameters(Parameters const* pars) const
{
  if(pars == NULL)
    pars = &(this->_pars);
  printf("MNR scale parameters:\n");
  printf("%f  %f  %f\n", pars->mf_A_c, pars->mf_B_c, pars->mf_C_c);
  printf("%f  %f  %f\n", pars->mr_A_c, pars->mr_B_c, pars->mr_C_c);
  printf("%f  %f  %f\n", pars->mf_A_b, pars->mf_B_b, pars->mf_C_b);
  printf("%f  %f  %f\n", pars->mr_A_b, pars->mr_B_b, pars->mr_C_b);
  printf("MNR masses:\n");
  printf("mc = %f  mb = %f\n", pars->mc, pars->mb);
  printf("MNR fragmentation parameters:\n");
  printf("fragpar_c = %f  fragpar_b = %f\n", pars->fragpar_c, pars->fragpar_b);
}
