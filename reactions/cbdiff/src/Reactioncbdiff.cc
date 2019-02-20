 
/*
   @file Reactioncbdiff.cc
   @date 2019-02-01
   @author  AddReaction.py
   Created by  AddReaction.py on 2019-02-01
*/

#include "Reactioncbdiff.h"
#include <TMath.h>
#include <TF1.h>

// the class factories
extern "C" Reactioncbdiff* create() {
  return new Reactioncbdiff();
}

void Reactioncbdiff::setDatasetParameters(int dataSetID, map<string,string> pars, map<string,double> dsPars)
{
  ReactionBaseHVQMNR::setDatasetParameters(dataSetID, pars, dsPars);
  //_debug = 1;
  Steering steer;
  if(checkParam("steer_q") || pars.find("steer_q") != pars.end())
    steer.q = GetParamIInPriority("steer_q", pars);
  else
    steer.q = true;
  if(checkParam("steer_a") || pars.find("steer_a") != pars.end())
    steer.a = GetParamIInPriority("steer_a", pars);
  else
    steer.a = true;
  steer.nf = GetParamIInPriority("steer_nf", pars);
  steer.ptmin = GetParamInPriority("steer_ptmin", pars);
  steer.ptmax = GetParamInPriority("steer_ptmax", pars);
  steer.npt = GetParamIInPriority("steer_npt", pars);
  steer.nptsm = GetParamIInPriority("steer_nptsm", pars);
  steer.ymin = GetParamInPriority("steer_ymin", pars);
  steer.ymax = GetParamInPriority("steer_ymax", pars);
  steer.ny = GetParamIInPriority("steer_ny", pars);
  steer.nsfnb = GetParamIInPriority("steer_nsfnb", pars);
  steer.nx3 = GetParamIInPriority("steer_nx3", pars);
  steer.nx4 = GetParamIInPriority("steer_nx4", pars);
  steer.nbz = GetParamIInPriority("steer_nbz", pars);
  steer.xmin = GetParamIInPriority("steer_xmin", pars);
  steer.xmax = GetParamIInPriority("steer_xmax", pars);
  steer.mf2min = GetParamIInPriority("steer_mf2min", pars);
  steer.mf2max = GetParamIInPriority("steer_mf2max", pars);

  std::shared_ptr<Parameters>& par = _mapPar[dataSetID];
  par = std::shared_ptr<Parameters>(new Parameters);
  // Order
  std::string order = GetParamSInPriority("Order", pars);
  MNR::MNRContribution contr = 11111; // NLO
  if(order == "LO")
    contr = 10111;
  else if(order != "NLO")
    hf_errlog(19020301, "F: order " + order + " not supported");
  const int ncontr = 1;
  MNR::MNRContribution** ptrContr = new MNR::MNRContribution*[ncontr];
  ptrContr[0] = new MNR::MNRContribution(contr);
  // HQ masses
  par->mc = GetParamInPriority("mq", pars);
  // scale parameters
  //GetMuPar('f', 'q', par->mf_A_c, par->mf_B_c, par->mf_C_c, pars);
  //GetMuPar('r', 'q', par->mr_A_c, par->mr_B_c, par->mr_C_c, pars);
  par->mf_A_c = GetParamInPriority("mf_A", pars);
  par->mf_B_c = GetParamInPriority("mf_B", pars);
  par->mr_A_c = GetParamInPriority("mr_A", pars);
  par->mr_B_c = GetParamInPriority("mr_B", pars);
  // fragmentation parameters
  par->fragpar_c = GetParamInPriority("FragPar", pars);
  PrintParameters(par.get());

  std::shared_ptr<MNR::MNR>& mnr = _mapMNR[dataSetID];
  mnr = std::shared_ptr<MNR::MNR>(new MNR::MNR(this));
  mnr->SetScaleCoef(par->mf_A_c, par->mf_B_c, par->mf_C_c, par->mr_A_c, par->mr_B_c, par->mr_C_c);

  std::shared_ptr<MNR::Grid>& grid = _mapGrid[dataSetID];
  grid = std::shared_ptr<MNR::Grid>(new MNR::Grid(ncontr, ptrContr));
  std::shared_ptr<MNR::Grid>& gridSm = _mapGridSm[dataSetID];
  gridSm = std::shared_ptr<MNR::Grid>(new MNR::Grid(ncontr, ptrContr));
  std::shared_ptr<MNR::Frag>& frag = _mapFrag[dataSetID];
  frag = std::shared_ptr<MNR::Frag>(new MNR::Frag);
  std::vector<TH2D*>& xsec = _mapXSec[dataSetID];

  mnr->SetDebug(_debug);
  DefaultInit(steer, par->mc, *mnr.get(), *frag.get(), *grid.get(), *gridSm.get());
  mnr->fC_sh = TMath::Power(stod(pars["energy"]), 2.0); // centre-of-mass energy squared
  mnr->CalcConstants();
  mnr->SetDebug(_debug);
  std::string finalState = pars["FinalState"];
  if(finalState == "parton")
  {
    frag->AddOut(NULL, par->mc);
  }
  else
  {
    frag->AddOut(MNR::Frag::GetFragFunction(0, finalState.c_str(), par->fragpar_c), MNR::Frag::GetHadronMass(finalState.c_str()));
    frag->GetFF(0)->SetParameter(1, par->fragpar_c);
  }
  _mapFF[dataSetID] = stod(pars["FragFrac"]);

  xsec.resize(1);
  for(size_t i = 0; i < xsec.size(); i++)
  {
    xsec[i] = new TH2D;
    // pT binning
    std::vector<double> binsPt;
    if(pars.find("pTn") != pars.end() && pars.find("pTmin") != pars.end() && pars.find("pTmax") != pars.end())
    {
      int nb = stoi(pars["pTn"]);
      binsPt.resize(nb + 1);
      double w = (stod(pars["pTmax"]) - stod(pars["pTmin"])) / nb;
      for(int b = 0; b < nb + 1; b++)
        binsPt[b] = stod(pars["pTmin"]) + w * b;
    }
    else if(pars.find("pT") != pars.end())
    {
      std::istringstream ss(pars["pT"]);
      std::string token;
      while(std::getline(ss, token, ','))
        binsPt.push_back(stod(token));
    }
    else
      hf_errlog(19021900, "F: no pT binning provided");
    // y binning
    std::vector<double> binsY;
    if(pars.find("yn") != pars.end() && pars.find("ymin") != pars.end() && pars.find("ymax") != pars.end())
    {
      int nb = stoi(pars["yn"]);
      binsY.resize(nb + 1);
      double w = (stod(pars["ymax"]) - stod(pars["ymin"])) / nb;
      for(int b = 0; b < nb + 1; b++)
        binsY[b] = stod(pars["ymin"]) + w * b;
    }
    else if(pars.find("y") != pars.end())
    {
      std::istringstream ss(pars["y"]);
      std::string token;
      while(std::getline(ss, token, ','))
        binsY.push_back(stod(token));
    }
    else
      hf_errlog(19021901, "F: no y binning provided");
    xsec[i]->SetBins(binsPt.size() - 1, &binsPt[0], binsY.size() - 1, &binsY[0]);
  }

  // test
  //mnr.get()->CalcXS(grid.get(), par->mc);
  //MNR::Grid::InterpolateGrid(gridSm.get(), gridSm.get(), par->mc);
  //frag->CalcCS(gridSm.get(), par->mc, xsec);
  //xsec[0]->Print("all");
}

// Initialize at the start of the computation
int Reactioncbdiff::initAtStart(const string &s)
{
  return 0;
}

void Reactioncbdiff::initAtIteration()
{
  ;
}

// Main function to compute results at an iteration
int Reactioncbdiff::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
  //printf("COMPUTE\n");
  std::shared_ptr<MNR::MNR> mnr(_mapMNR[dataSetID]);
  std::shared_ptr<Parameters> par(_mapPar[dataSetID]);
  std::shared_ptr<MNR::Grid> grid(_mapGrid[dataSetID]);
  std::shared_ptr<MNR::Grid> gridSm(_mapGridSm[dataSetID]);
  std::shared_ptr<MNR::Frag> frag(_mapFrag[dataSetID]);
  std::vector<TH2D*> xsec(_mapXSec[dataSetID]);

  mnr->CalcXS(grid.get(), par->mc);
  MNR::Grid::InterpolateGrid(grid.get(), gridSm.get(), par->mc);
  frag->CalcCS(gridSm.get(), par->mc, xsec);
  //xsec[0]->Print("all");

  DataSet& ds = _dataSets[dataSetID];
  if (ds.BinsYMin == NULL || ds.BinsYMax == NULL || ds.BinsPtMin == NULL || ds.BinsPtMax == NULL )
  {
    // fill results array with cross sections bin by bin
    int nbx = xsec[0]->GetNbinsX();
    for(size_t i = 0; i < val.size(); i++)
      val[i] = xsec[0]->GetBinContent((i % nbx) + 1, (i / nbx) + 1) * _mapFF[dataSetID];
  }
  else
  {
    // fill results array with cross sections by matching bins
    for(unsigned int i = 0; i < ds.BinsYMin->size(); i++)
      val[i] = FindXSecPtYBin(xsec[0], (*ds.BinsYMin)[i], (*ds.BinsYMax)[i], (*ds.BinsPtMin)[i], (*ds.BinsPtMax)[i], false, false) * _mapFF[dataSetID];
  }

  return 0;
}

