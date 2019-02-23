 
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
  steer.xmin = GetParamInPriority("steer_xmin", pars);
  steer.xmax = GetParamInPriority("steer_xmax", pars);
  steer.mf2min = GetParamInPriority("steer_mf2min", pars);
  steer.mf2max = GetParamInPriority("steer_mf2max", pars);

  std::shared_ptr<Parameters>& par = _mapPar[dataSetID];
  par = std::shared_ptr<Parameters>(new Parameters);
  // Order
  std::string order = GetParamSInPriority("Order", pars);
  MNR::MNRContribution contrNLO = 11111; // NLO
  MNR::MNRContribution contrLO = 10111; // NLO
  MNR::MNRContribution contr = contrNLO;
  if(order == "LO")
    contr = contrLO;
  else if(order != "NLO")
    hf_errlog(19020301, "F: order " + order + " not supported");
  const int ncontr = 1;
  MNR::MNRContribution** ptrContr = new MNR::MNRContribution*[ncontr];
  ptrContr[0] = new MNR::MNRContribution(contr);
  MNR::MNRContribution** ptrContrLO = new MNR::MNRContribution*[ncontr];
  ptrContrLO[0] = new MNR::MNRContribution(contrLO);
  // HQ masses
  par->mc = GetParamInPriority("mq", pars);
  _mapMassDiff[dataSetID] = 0.001; // for MSbar mass transformation; 1 MeV should work for c, b and t
  //_mapMassDiff[dataSetID] = 0.150; // for MSbar mass transformation; 1 MeV should work for c, b and t
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
  // pole or MSbar mass
  _mapMSbarMass[dataSetID] = 0;
  if(checkParamInPriority("MS_MASS", pars))
    _mapMSbarMass[dataSetID] = GetParamInPriority("MS_MASS", pars);

  std::shared_ptr<MNR::MNR>& mnr = _mapMNR[dataSetID];
  mnr = std::shared_ptr<MNR::MNR>(new MNR::MNR(this));
  mnr->SetScaleCoef(par->mf_A_c, par->mf_B_c, par->mf_C_c, par->mr_A_c, par->mr_B_c, par->mr_C_c);

  // main grid (LO or NLO)
  std::shared_ptr<MNR::Grid>& grid = _mapGrid[dataSetID];
  grid = std::shared_ptr<MNR::Grid>(new MNR::Grid(ncontr, ptrContr));
  // LO mass up grid (for MSbar mass)
  std::shared_ptr<MNR::Grid>& gridLOMassUp = _mapGridLOMassUp[dataSetID];
  gridLOMassUp = std::shared_ptr<MNR::Grid>(new MNR::Grid(ncontr, ptrContrLO));
  std::shared_ptr<MNR::Grid>& gridSmLOMassUp = _mapGridSmLOMassUp[dataSetID];
  gridSmLOMassUp = std::shared_ptr<MNR::Grid>(new MNR::Grid(ncontr, ptrContrLO));
  // LO mass down grid (for MSbar mass)
  std::shared_ptr<MNR::Grid>& gridLOMassDown = _mapGridLOMassDown[dataSetID];
  gridLOMassDown = std::shared_ptr<MNR::Grid>(new MNR::Grid(ncontr, ptrContrLO));
  std::shared_ptr<MNR::Grid>& gridSmLOMassDown = _mapGridSmLOMassDown[dataSetID];
  gridSmLOMassDown = std::shared_ptr<MNR::Grid>(new MNR::Grid(ncontr, ptrContrLO));
  // final smoothed grid
  std::shared_ptr<MNR::Grid>& gridSm = _mapGridSm[dataSetID];
  gridSm = std::shared_ptr<MNR::Grid>(new MNR::Grid(ncontr, ptrContr));

  std::shared_ptr<MNR::Frag>& frag = _mapFrag[dataSetID];
  frag = std::shared_ptr<MNR::Frag>(new MNR::Frag);
  std::vector<TH2D*>& xsec = _mapXSec[dataSetID];

  mnr->SetDebug(_debug);
  DefaultInitMNR(steer, par->mc, *mnr.get());
  DefaultInitGrid(steer, par->mc, steer.npt, *grid.get());
  DefaultInitGrid(steer, par->mc, steer.npt, *gridLOMassUp.get());
  DefaultInitGrid(steer, par->mc, steer.npt, *gridLOMassDown.get());
  //DefaultInitGrid(steer, par->mc + _mapMassDiff[dataSetID], steer.npt, *gridLOMassUp.get());
  //DefaultInitGrid(steer, par->mc - _mapMassDiff[dataSetID], steer.npt, *gridLOMassDown.get());
  DefaultInitGrid(steer, par->mc, steer.nptsm, *gridSm.get());
  DefaultInitGrid(steer, par->mc, steer.nptsm, *gridSmLOMassUp.get());
  DefaultInitGrid(steer, par->mc, steer.nptsm, *gridSmLOMassDown.get());
  DefaultInitFrag(steer, *frag.get());
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
  //MNR::Grid::InterpolateGrid(grid.get(), gridSm.get(), par->mc);

  // tarnsformation to MSbar mass scheme
  if(_mapMSbarMass[dataSetID])
  {
    // store original scale B parameters which need to be modify for changed mass
    const double mfB = mnr->fMf_B;
    const double mrB = mnr->fMr_B;

    // LO mass up variation
    double massU = par->mc + _mapMassDiff[dataSetID];
    mnr->fMf_B = mfB * pow(par->mc / massU, 2.0);
    mnr->fMr_B = mrB * pow(par->mc / massU, 2.0);
    std::shared_ptr<MNR::Grid> gridLOMassU(_mapGridLOMassUp[dataSetID]);
    mnr->CalcXS(gridLOMassU.get(), massU);
    //std::shared_ptr<MNR::Grid> gridSmLOMassU(_mapGridSmLOMassUp[dataSetID]);
    //MNR::Grid::InterpolateGrid(gridLOMassU.get(), gridSmLOMassU.get(), massU);
    //MNR::Grid::InterpolateGrid(gridLOMassU.get(), gridSmLOMassU.get(), par->mc);

    // LO mass down variation
    double massD = par->mc - _mapMassDiff[dataSetID];
    mnr->fMf_B = mfB * pow(par->mc / massD, 2.0);
    mnr->fMr_B = mrB * pow(par->mc / massD, 2.0);
    std::shared_ptr<MNR::Grid> gridLOMassD(_mapGridLOMassDown[dataSetID]);
    mnr->CalcXS(gridLOMassD.get(), massD);
    //std::shared_ptr<MNR::Grid> gridSmLOMassD(_mapGridSmLOMassDown[dataSetID]);
    //MNR::Grid::InterpolateGrid(gridLOMassD.get(), gridSmLOMassD.get(), massD);
    //MNR::Grid::InterpolateGrid(gridLOMassD.get(), gridSmLOMassD.get(), par->mc);

    // scheme transformation
    //MNR::Grid::TransformGridToMSbarMassScheme(grid.get(), gridLOMassU.get(), gridLOMassD.get(), par->mc, _mapMassDiff[dataSetID]);
    //MNR::Grid::TransformGridToMSbarMassScheme(gridSm.get(), gridSmLOMassU.get(), gridSmLOMassD.get(), par->mc, _mapMassDiff[dataSetID]);

    // restore original scales
    mnr->fMf_B = mfB;
    mnr->fMr_B = mrB;

    MNR::Grid::InterpolateGrid(grid.get(), gridSm.get(), par->mc, gridLOMassU.get(), massU, gridLOMassD.get(), massD);
    //MNR::Grid::InterpolateGrid(grid.get(), gridSm.get(), par->mc);
  }
  else
    MNR::Grid::InterpolateGrid(grid.get(), gridSm.get(), par->mc);

  //MNR::Grid::InterpolateGrid(grid.get(), gridSm.get(), par->mc);
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

