
/*
   @file Reactioncbdiff.cc
   @date 2019-02-01
   @author  AddReaction.py
   Created by  AddReaction.py on 2019-02-01
*/

#include "Reactioncbdiff.h"
#include <TMath.h>
#include <TF1.h>
#include "xfitter_cpp_base.h"

// the class factories
extern "C" Reactioncbdiff* create() {
  return new Reactioncbdiff();
}

void Reactioncbdiff::initTerm(TermData *td)
{
  ReactionBaseHVQMNR::initTerm(td);
  unsigned dataSetID = td->id;

  //_debug = 1;
  Steering steer;
  steer.q = td->hasParam("steer_q") ? td->getParamI("steer_q") : true;
  steer.a = td->hasParam("steer_a") ? td->getParamI("steer_a") : true;
  steer.nf = td->getParamI("steer_nf");
  steer.ptmin = *td->getParamD("steer_ptmin");
  steer.ptmax = *td->getParamD("steer_ptmax");
  steer.npt = td->getParamI("steer_npt");
  steer.nptsm = td->getParamI("steer_nptsm");
  steer.ymin = *td->getParamD("steer_ymin");
  steer.ymax = *td->getParamD("steer_ymax");
  steer.ny = td->getParamI("steer_ny");
  steer.nsfnb = td->getParamI("steer_nsfnb");
  steer.nx3 = td->getParamI("steer_nx3");
  steer.nx4 = td->getParamI("steer_nx4");
  steer.nbz = td->getParamI("steer_nbz");
  steer.xmin = *td->getParamD("steer_xmin");
  steer.xmax = *td->getParamD("steer_xmax");
  steer.mf2min = *td->getParamD("steer_mf2min");
  steer.mf2max = *td->getParamD("steer_mf2max");

  // precision: 1.0 is default
  _mapPrecision[dataSetID] = 1.0;
  if(td->hasParam("precision"))
  {
    _mapPrecision[dataSetID] = *td->getParamD("precision");
    if(_mapPrecision[dataSetID] != 1.0)
    {
      printf("Using precision factor %f\n", _mapPrecision[dataSetID]);
      steer.npt *= _mapPrecision[dataSetID];
      steer.nptsm *= _mapPrecision[dataSetID];
      steer.ny *= _mapPrecision[dataSetID];
      steer.nsfnb *= _mapPrecision[dataSetID];
      steer.nx3 *= _mapPrecision[dataSetID];
      steer.nx4 *= _mapPrecision[dataSetID];
      steer.nbz *= _mapPrecision[dataSetID];
    }
  }

  std::shared_ptr<Parameters>& par = _mapPar[dataSetID];
  par = std::shared_ptr<Parameters>(new Parameters);
  // Flavour
  par->flav = td->hasParam("flav") ? td->getParamS("flav").c_str()[0] : 'c';
  // Order
  const string order = td->getParamS("Order");
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
  if(td->hasParam("mq"))
    par->mc = *td->getParamD("mq");
  else
  {
    par->flagMassIsGlobal = true;
    if(par->flav == 'c')
      par->mc = *td->getParamD("mch");
    else if(par->flav == 'b')
      par->mc = *td->getParamD("mbt");
    else if(par->flav == 't')
      par->mc = *td->getParamD("mtp");
  }
  _mapMassDiff[dataSetID] = 0.001; // for MSbar mass transformation; 1 MeV should work for c, b and t
  //_mapMassDiff[dataSetID] = 0.150;
  // scale parameters
  //GetMuPar('f', 'q', par->mf_A_c, par->mf_B_c, par->mf_C_c, pars);
  //GetMuPar('r', 'q', par->mr_A_c, par->mr_B_c, par->mr_C_c, pars);
  par->mf_A_c = *td->getParamD("mf_A");
  par->mf_B_c = *td->getParamD("mf_B");
  par->mf_C_c = 0.0;
  if(td->hasParam("mf_C"))
    par->mf_C_c = *td->getParamD("mf_C");
  par->mr_A_c = *td->getParamD("mr_A");
  par->mr_B_c = *td->getParamD("mr_B");
  par->mr_C_c = 0.0;
  if(td->hasParam("mr_C"))
    par->mr_C_c = *td->getParamD("mr_C");
  // fragmentation parameters
  par->fragpar_c = *td->getParamD("FragPar");
  PrintParameters(par.get());
  // pole, MSbar or MSR mass (0 pole, 1 MSbar, 2 MSR, 10 MSbar with log term in d1, 11 MSbar with varied mu_m - mq_mu0 must be secified)
  _mapMassScheme[dataSetID] = td->hasParam("MS_MASS") ? td->getParamI("MS_MASS") : 0;
  printf("MNR: order = %s MS_MASS = %d\n", order.c_str(), _mapMassScheme[dataSetID]);
  // divide or not by bin width
  par->flagDivideBinWidth = td->hasParam("dividebw") ? td->getParamI("dividebw") : false;
  // debug mode
  if(td->hasParam("debug"))
  par->debug = td->hasParam("dividebw") ? td->getParamI("debug") : 0;

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
  DefaultInitGrid(steer, par->mc, steer.nptsm, *gridSm.get());
  DefaultInitGrid(steer, par->mc, steer.nptsm, *gridSmLOMassUp.get());
  DefaultInitGrid(steer, par->mc, steer.nptsm, *gridSmLOMassDown.get());
  DefaultInitFrag(steer, *frag.get());
  mnr->fC_sh = TMath::Power(*td->getParamD("energy"), 2.0); // centre-of-mass energy squared
  mnr->CalcConstants();
  std::string finalState = td->getParamS("FinalState");
  if(finalState == "parton")
  {
    frag->AddOut(NULL, par->mc);
    _mapHadronMass[dataSetID] = 0.0; // hadron mass is irrelevant
  }
  else
  {
    int fragType = 0; // Kartvelishvili
    // read hadron mass (if <0 then PDG values are used)
    _mapHadronMass[dataSetID] = td->hasParam("hadronMass") ? *td->getParamD("hadronMass") : -1.0;
    if(_mapHadronMass[dataSetID] < 0.0)
      _mapHadronMass[dataSetID] = MNR::Frag::GetHadronMass(finalState.c_str());
    //frag->AddOut(MNR::Frag::GetFragFunction(fragType, finalState.c_str(), par->fragpar_c), MNR::Frag::GetHadronMass(finalState.c_str()));
    printf("MNRFrag: using Kartvelishvili function with par = %.1f and hadronMass = %.3f\n", par->fragpar_c, _mapHadronMass[dataSetID]);
    frag->AddOut(MNR::Frag::GetFragFunction(fragType, finalState.c_str(), par->fragpar_c), _mapHadronMass[dataSetID]);
    frag->GetFF(0)->SetParameter(1, par->fragpar_c);
  }
  _mapFF[dataSetID] = *td->getParamD("FragFrac");

  xsec.resize(1);
  for(size_t i = 0; i < xsec.size(); i++)
  {
    xsec[i] = new TH2D;
    // pT binning
    std::vector<double> binsPt;
    if(td->hasParam("pTn") && td->hasParam("pTmin") && td->hasParam("pTmax"))
    {
      int nb = td->getParamI("pTn");
      binsPt.resize(nb + 1);
      double w = (*td->getParamD("pTmax") - *td->getParamD("pTmin")) / nb;
      for(int b = 0; b < nb + 1; b++)
        binsPt[b] = *td->getParamD("pTmin") + w * b;
    }
    else if(td->hasParam("pT"))
    {
      std::istringstream ss(td->getParamS("pT"));
      std::string token;
      while(std::getline(ss, token, ','))
        binsPt.push_back(stod(token));
    }
    else
      hf_errlog(19021900, "F: no pT binning provided");
    // y binning
    std::vector<double> binsY;
    if(td->hasParam("yn") && td->hasParam("ymin") && td->hasParam("ymax"))
    {
      int nb = td->getParamI("yn");
      binsY.resize(nb + 1);
      double w = (*td->getParamD("ymax") - *td->getParamD("ymin")) / nb;
      for(int b = 0; b < nb + 1; b++)
        binsY[b] = *td->getParamD("ymin") + w * b;
    }
    else if(td->hasParam("y"))
    {
      std::istringstream ss(td->getParamS("y"));
      std::string token;
      while(std::getline(ss, token, ','))
        binsY.push_back(stod(token));
    }
    else
      hf_errlog(19021901, "F: no y binning provided");
    xsec[i]->SetBins(binsPt.size() - 1, &binsPt[0], binsY.size() - 1, &binsY[0]);
  }
  _mapN[dataSetID] = 1;
  if(td->hasParam("N"))
    _mapN[dataSetID] = td->getParamI("N");
}

// Main function to compute results at an iteration
void Reactioncbdiff::compute(TermData *td, valarray<double> &val, map<string, valarray<double>> &errors)
{
  td->actualizeWrappers();
  int dataSetID = td->id;
  std::shared_ptr<MNR::MNR> mnr(_mapMNR[dataSetID]);
  std::shared_ptr<Parameters> par(_mapPar[dataSetID]);
  std::shared_ptr<MNR::Grid> grid(_mapGrid[dataSetID]);
  std::shared_ptr<MNR::Grid> gridSm(_mapGridSm[dataSetID]);
  std::shared_ptr<MNR::Frag> frag(_mapFrag[dataSetID]);
  std::vector<TH2D*> xsec(_mapXSec[dataSetID]);

  // update mass which may be fitted (scales or frag. par. cannot be fitted with this implementation)
  if(par->flagMassIsGlobal)
  {
    if(par->flav == 'c')
      par->mc = *td->getParamD("mch");
    else if(par->flav == 'b')
      par->mc = *td->getParamD("mbt");
    else if(par->flav == 't')
      par->mc = *td->getParamD("mtp");
  }

  // protection against not positive or nan mass
  if(par->mc < 0.0 || par->mc != par->mc)
    par->mc = 1000.0; // probably there are no heavy quarks for which mass of 1000 GeV would be reasonable

  double mc_mu0 = -1.0;
  double mc_mu0_var = 0;
  if(_mapMassScheme[dataSetID] == 11) 
  {
    mc_mu0 = *td->getParamD("mq_mu0");
    mc_mu0_var = *td->getParamD("mq_mu0_var");
  }
  mnr->CalcXS(grid.get(), par->mc, mc_mu0);

  // tarnsformation to MSbar mass scheme
  if(_mapMassScheme[dataSetID])
  {
    // store original scale B parameters which need to be modified for changed mass
    const double mfB = mnr->fMf_B;
    const double mrB = mnr->fMr_B;

    // LO mass up variation
    double massU = par->mc + _mapMassDiff[dataSetID];
    mnr->fMf_B = mfB * pow(par->mc / massU, 2.0);
    mnr->fMr_B = mrB * pow(par->mc / massU, 2.0);
    std::shared_ptr<MNR::Grid> gridLOMassU(_mapGridLOMassUp[dataSetID]);
    mnr->CalcXS(gridLOMassU.get(), massU, mc_mu0);

    // LO mass down variation
    double massD = par->mc - _mapMassDiff[dataSetID];
    mnr->fMf_B = mfB * pow(par->mc / massD, 2.0);
    mnr->fMr_B = mrB * pow(par->mc / massD, 2.0);
    std::shared_ptr<MNR::Grid> gridLOMassD(_mapGridLOMassDown[dataSetID]);
    mnr->CalcXS(gridLOMassD.get(), massD, mc_mu0);

    // restore original scales
    mnr->fMf_B = mfB;
    mnr->fMr_B = mrB;

    if(_mapMassScheme[dataSetID] == 1)
    {
      int flagMSbarTransformation = 0; // d1=4/3 (no ln)
      MNR::Grid::InterpolateGrid(grid.get(), gridSm.get(), par->mc, gridLOMassU.get(), massU, gridLOMassD.get(), massD, flagMSbarTransformation);
    }
    else if(_mapMassScheme[dataSetID] == 10)
    {
      int flagMSbarTransformation = 1; // d1=4/3 + ln
      MNR::Grid::InterpolateGrid(grid.get(), gridSm.get(), par->mc, gridLOMassU.get(), massU, gridLOMassD.get(), massD, flagMSbarTransformation);
    }
    else if(_mapMassScheme[dataSetID] == 2)
    {
      double R = 3.0;
      if(par->flav == 'c')
        R = 1.0;
      int nl = mnr->GetNl();
      MNR::Grid::InterpolateGrid(grid.get(), gridSm.get(), par->mc, gridLOMassU.get(), massU, gridLOMassD.get(), massD, 2, &R, &nl);
    }
    else if(_mapMassScheme[dataSetID] == 11)
    {
      // MSbar with mu_m variation as used in 2009.07763
      MNR::Grid::InterpolateGrid(grid.get(), gridSm.get(), par->mc, gridLOMassU.get(), massU, gridLOMassD.get(), massD, 11, 0, 0, mc_mu0*mc_mu0_var);
    }
  }
  else
    MNR::Grid::InterpolateGrid(grid.get(), gridSm.get(), par->mc);

  frag->CalcCS(gridSm.get(), par->mc, xsec);
  if(par->debug)
    xsec[0]->Print("all");

  DataSet& ds = _dataSets[dataSetID];
  val.resize(xsec[0]->GetNbinsX() * xsec[0]->GetNbinsY() * _mapN[dataSetID]);
  if (ds.BinsYMin == NULL || ds.BinsYMax == NULL || ds.BinsPtMin == NULL || ds.BinsPtMax == NULL )
  {
    // fill results array with cross sections bin by bin, optionally repeat N times
    for(int in = 0; in < _mapN[dataSetID]; in++)
    {
      for(int ix = 1; ix <= xsec[0]->GetNbinsX(); ix++)
      {
        for(int iy = 1; iy <= xsec[0]->GetNbinsY(); iy++)
        {
          //int ival = (ix - 1) * xsec[0]->GetNbinsY() + iy - 1 + in * xsec[0]->GetNbinsX() * xsec[0]->GetNbinsY();
          int ival = (iy - 1) * xsec[0]->GetNbinsX() + ix - 1 + in * xsec[0]->GetNbinsX() * xsec[0]->GetNbinsY();
          val[ival] = xsec[0]->GetBinContent(ix, iy) * _mapFF[dataSetID];
          if(par->flagDivideBinWidth)
          {
            val[ival] /= xsec[0]->GetXaxis()->GetBinWidth(ix);
            val[ival] /= xsec[0]->GetYaxis()->GetBinWidth(iy);
          }
        }
      }
    }
  }
  else
  {
    // fill results array with cross sections by matching bins
    // make sure number of bins is consistent
    if(ds.BinsYMin->size() > (size_t)(xsec[0]->GetNbinsX() * xsec[0]->GetNbinsY()))
      hf_errlog(19021900, "F: inconsistent number of bins");
    for(unsigned int i = 0; i < ds.BinsYMin->size(); i++)
      val[i] = FindXSecPtYBin(xsec[0], (*ds.BinsYMin)[i], (*ds.BinsYMax)[i], (*ds.BinsPtMin)[i], (*ds.BinsPtMax)[i], par->flagDivideBinWidth, par->flagDivideBinWidth) * _mapFF[dataSetID];
  }
  if(par->debug)
  {
    for(size_t i = 0; i < val.size(); i++)
      printf("val[%lu] = %f\n", i, val[i]);
  }
}

