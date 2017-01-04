#include "xfitter_cpp.h"
#include <TMath.h>
 
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
  std::string str = "F: can not create ReactionBaseHVQMNR instance: you should implement calculation in a derived class";
  hf_errlog_(16123000, str.c_str(), str.length());
  return NULL;
}


// pass to MNR pointer to instance inherited from ReactionTheory to allow access to alphas and PDF routines
ReactionBaseHVQMNR::ReactionBaseHVQMNR() : _mnr(MNR::MNR(this))
{
  //printf("OZ ReactionBaseHVQMNR::ReactionBaseHVQMNR()\n");
  // set initialisation status flag
  _isInitAtStart = false;
  // set debugging flag
  _debug = steering_.ldebug_;
  //_debug = 1;
}


ReactionBaseHVQMNR::~ReactionBaseHVQMNR()
{
  //printf("OZ ReactionBaseHVQMNR::~ReactionBaseHVQMNR()\n");
  for(unsigned int i = 0; i < _hCalculatedXSec.size(); i++)
    delete _hCalculatedXSec[i];
}


void ReactionBaseHVQMNR::setDatasetParamters(int dataSetID, map<string,string> pars)
{
  // add new dataset
  DataSet dataSet;
  std::pair<std::map<int, DataSet>::iterator, bool> ret = _dataSets.insert(std::pair<int, DataSet>(dataSetID, dataSet));
  // check if dataset with provided ID already exists
  if(!ret.second)
  {
    char buffer[256];
    sprintf(buffer, "F: dataset with id = %d already exists", dataSetID);
    error(16123001, buffer);
  }

  // set parameters for new dataset
  std::string str = pars.begin()->first;
  DataSet& ds = ret.first->second;
  // "FinalState=" must be provided
  if(readFromTermInfo(str, "FinalState=", ds.FinalState))
    error(16123002, "F: TermInfo must contain FinalState=...");

  // optional "NormY="
  if(readFromTermInfo(str, "NormY=", ds.NormY))
    ds.NormY = 0; // default value is unnormalised absolute cross section

  // optional "FragFrac=" (must be provided for absolute cross section)
  if(readFromTermInfo(str, "FragFrac=", ds.FragFraction))
  {
    if(ds.NormY == 0)
      error(16123003, "F: for absolute cross section TermInfo must contain FragFrac=...");
  }
  else
  {
    if(ds.NormY != 0)
      printf("Warning: FragFrac=%f will be ignored for normalised cross sections\n", ds.FragFraction);
  }
  
  // set binning
  ds.BinsYMin  = GetBinValues(dataSetID, "ymin");
  ds.BinsYMax  = GetBinValues(dataSetID, "ymax");
  ds.BinsPtMin = GetBinValues(dataSetID, "pTmin");
  ds.BinsPtMax = GetBinValues(dataSetID, "pTmax");
  if (ds.BinsYMin == NULL || ds.BinsYMax == NULL || ds.BinsPtMin == NULL || ds.BinsPtMax == NULL ) 
    error(16123004, "F: No bins ymin or ymax or ptmin or ptmax");
  // set reference y bins if needed
  if(ds.NormY == 1)
  {
    ds.BinsYMinRef = GetBinValues(dataSetID, "yminREF");
    ds.BinsYMaxRef = GetBinValues(dataSetID, "ymaxREF");
    if(ds.BinsYMinRef == NULL || ds.BinsYMaxRef == NULL)
      error(16123005, "F: No bins yminREF or ymaxREF for normalised cross section");
  }
  
  if(_debug)
    printf("Added dataset: FinalState = %s  NormY = %d  FragFraction = %f\n", ds.FinalState.c_str(), ds.NormY, ds.FragFraction);
}


// ********************************
// ***** utility routines *********
// ********************************

// error treatment
void ReactionBaseHVQMNR::error(const int id, const std::string& str)
{
  hf_errlog_(id, str.c_str(), str.length());
}

// check equality of float numbers with tolerance
bool ReactionBaseHVQMNR::IsEqual(const double val1, const double val2, const double eps/* = 1e-6*/)
{
  return (TMath::Abs(val1 - val2) < eps || TMath::Abs((val1 - val2) / val1) < eps);
}

// initialise calculation with default parameters
void ReactionBaseHVQMNR::DefaultInit(const Steering& steer, const double mq, MNR::MNR& mnr, MNR::Frag& frag, MNR::Grid& grid, MNR::Grid& grid_smoothed)
{
  // MNR (parton level cross sections)
  mnr.bFS_Q = true;
  mnr.bFS_A = true;
  // x3 and x4 binning
  mnr.fBn_x3 = steer.nx3;
  mnr.fBn_x4 = steer.nx4;
  mnr.fSF_nb = steer.nsfnb;
  mnr.CalcBinning();
  // Number of flavours
  mnr.fC_nl = 3;
  // Parton level pT-y grids
  grid.SetL(steer.npt, steer.ptmin, steer.ptmax, mq);
  grid.SetY(steer.ny, steer.ymin, steer.ymax);
  grid.SetW(1);
  grid_smoothed.SetL(steer.nptsm, steer.ptmin, steer.ptmax, mq);
  grid_smoothed.SetY(steer.ny, steer.ymin, steer.ymax);
  grid_smoothed.SetW(1);
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
  char str[256];
  sprintf(str, "F: ERROR in FindXSPtY(): bin not found y: %f %f, pt: %f %f\n", ymin, ymax, ptmin, ptmax);
  error(16123006, str);
  return -1.0;
}

// check if appropriate heavy-flavour scheme is used
void ReactionBaseHVQMNR::CheckHFScheme()
{
  // check HF scheme
  if(steering_.hfscheme_ != 3 && steering_.hfscheme_ != 4)
  {
    char str[256];
    sprintf(str, "S: calculation does not support HFSCHEME = %d (only 3, 4 supported)", steering_.hfscheme_);
    error(16123007, str);
  }
}

// read parameters for perturbative scales from MINUIT extra parameters
void ReactionBaseHVQMNR::GetMuPar(const char mu, const char q, double& A, double& B, double& C)
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
  if(checkParam(baseParameterName + "_AB"))
    A = B = GetParam(baseParameterName + "_AB");
  else
  {
    if(checkParam(baseParameterName + "_A") && checkParam(baseParameterName + "_B"))
    {
      A = GetParam(baseParameterName + "_A");
      B = GetParam(baseParameterName + "_B");
    }
    else
    {
      if(checkParam(baseParameterName + "_AB_" + std::string(1, q)))
        A = B = GetParam(baseParameterName + "_AB_" + std::string(1, q));
      else
      {
        if(checkParam(baseParameterName + "_A_" + std::string(1, q)))
          A = GetParam(baseParameterName + "_A_" + std::string(1, q));
        else
          A = defA;
        if(checkParam(baseParameterName + "_B_" + std::string(1, q)))
          B = GetParam(baseParameterName + "_B_" + std::string(1, q));
        else
          B = defB;
      }
    }
  }
  // C parameter
  if(checkParam(baseParameterName + "_C"))
    C = GetParam(baseParameterName + "_C");
  else
  {
    if(checkParam(baseParameterName + "_C_" + std::string(1, q)))
      C = GetParam(baseParameterName + "_C_" + std::string(1, q));
    else
      C = defC;
  }
}


// read fragmentation parameter from MINUIT extra parameters
double ReactionBaseHVQMNR::GetFragPar(const char q)
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
  double parvalue;
  char parname[16];
  sprintf(parname, "MNRfrag_%c", q);
  if(!checkParam(parname))
  {
    // parameter not in ExtraParamMinuit -> using default value
    if(q == 'c') 
      parvalue = defFFc;
    else if(q == 'b')
      parvalue = defFFb;
    else
      printf("Warning in GetFragPar(): no default value for q = %c\n", q);
    parvalue = NAN;
  }
  else
    parvalue = GetParam(parname);
    
  // TODO check below
  /*      ! parameter in ExtraParamMinuit, but not in MINUIT: this happens, if we are not in 'Fit' mode -> using default value
        if(st.lt.0) then
          if(q.eq.'c') then
            FFpar=defFFc
          else if(q.eq.'b') then
            FFpar=defFFb
          else
            write(*,*)'Warning in GetFPar(): no default value for q = ',q
            call makenan(FFpar)
          endif
        endif
      endif
      end*/
  return parvalue;
}

// read and update theory parameters
void ReactionBaseHVQMNR::UpdateParameters()
{
  // heavy-quark masses
  _pars.mc = fermion_masses_.mch_;
  _pars.mb = fermion_masses_.mbt_;
  // scale parameters
  GetMuPar('f', 'c', _pars.mf_A_c, _pars.mf_B_c, _pars.mf_C_c);
  GetMuPar('r', 'c', _pars.mr_A_c, _pars.mr_B_c, _pars.mr_C_c);
  GetMuPar('f', 'b', _pars.mf_A_b, _pars.mf_B_b, _pars.mf_C_b);
  GetMuPar('r', 'b', _pars.mr_A_b, _pars.mr_B_b, _pars.mr_C_b);
  // fragmentation parameters
  _pars.fragpar_c = GetFragPar('c');
  _pars.fragpar_b = GetFragPar('b');
}

// print theory parameters
void ReactionBaseHVQMNR::PrintParameters() const
{
  printf("MNR scale parameters:\n");
  printf("%f  %f  %f\n", _pars.mf_A_c, _pars.mf_B_c, _pars.mf_C_c);
  printf("%f  %f  %f\n", _pars.mr_A_c, _pars.mr_B_c, _pars.mr_C_c);
  printf("%f  %f  %f\n", _pars.mf_A_b, _pars.mf_B_b, _pars.mf_C_b);
  printf("%f  %f  %f\n", _pars.mr_A_b, _pars.mr_B_b, _pars.mr_C_b);
  printf("MNR masses:\n");
  printf("mc = %f  mb = %f\n", _pars.mc, _pars.mb);
  printf("MNR fragmentation parameters:\n");
  printf("fragpar_c = %f  fragpar_b = %f\n", _pars.fragpar_c, _pars.fragpar_b);
}

// read string value for provided key
int ReactionBaseHVQMNR::readFromTermInfo(const std::string& str, const std::string& key, std::string& value)
{
  //printf("read str = %s\n", str.c_str());
  std::size_t pos1 = str.find(key);
  // key not found
  if(pos1 == std::string::npos)
    return 1;
  pos1 += key.length(); // skip key length
  std::size_t pos2 = str.find(":", pos1);
  if(pos2 != std::string::npos)
    // key=value not in the end
    value = std::string(str, pos1, pos2 - pos1);
  else
    // key=value in the end
    value = std::string(str, pos1);
  return 0;
}

// read int value for provided key
int ReactionBaseHVQMNR::readFromTermInfo(const std::string& str, const std::string& key, int& value)
{
  std::string valueString;
  readFromTermInfo(str, key, valueString);
  if(valueString.length() == 0)
    return 1;
  value = atoi(valueString.c_str());
  return 0;
}

// read float value for provided key
int ReactionBaseHVQMNR::readFromTermInfo(const std::string& str, const std::string& key, float& value)
{
  std::string valueString;
  readFromTermInfo(str, key, valueString);
  if(valueString.length() == 0)
    return 1;
  value = atof(valueString.c_str());
  return 0;
}

// read double value for provided key
int ReactionBaseHVQMNR::readFromTermInfo(const std::string& str, const std::string& key, double& value)
{
  float value_float = 0.0;
  int status = readFromTermInfo(str, key, value_float);
  value = value_float;
  return status;
}
