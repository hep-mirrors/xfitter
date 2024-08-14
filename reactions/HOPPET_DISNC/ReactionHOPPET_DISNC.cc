/*
   @file Reaction_DISNC_Hoppet.cc
   @date 2017-04-08
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-04-08
*/

#include "ReactionHOPPET_DISNC.h"
#include <iostream>
#include <cstdio>
#include "ReactionBaseDISNC.h"
#include "hoppet_v1.h" // Include the HOPPET header
//#include <xfitter_cpp.h>
#include "xfitter_pars.h"
//#include "xfitter_cpp_base.h"
#include <BaseEvolution.h>
#include <hf_errlog.h>



using namespace hoppetv1;
// the class factories
extern "C" ReactionHOPPET_DISNC *create()
{
  return new ReactionHOPPET_DISNC();
}

// Initialize at the start of the computation
//void Reaction_DISNC_Hoppet::atStart()
//{
//}

// Main function to compute results at an iteration
void ReactionHOPPET_DISNC::F2(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal)
{
	const double mc = 1.414213563;
  const double mb = 4.5;
  const double mt = 175.0;
  hoppetSetPoleMassVFN(mc,mb,mt);
    
  const bool param_coefs = true;
  int	 order_max = 4;
	const double xmuR = 1.;
  const double xmuF = 1.; 
  hoppetInitStrFct(order_max, param_coefs, xmuR, xmuF);
    
  const auto _convfac = *XFITTER_PARS::getParamD("convFac");
  const auto _alphaem = *XFITTER_PARS::getParamD("alphaem");
  const auto MZ = *XFITTER_PARS::getParamD("Mz");
  const auto MW = *XFITTER_PARS::getParamD("Mw");
  //const auto _sin2thetaW = *XFITTER_PARS::getParamD("sin2thW");
  
  const double MW2 = MW*MW; 
  const double MZ2 = MZ*MZ;

  // Construct structure functions
  const double s2tw = 1 - MW2 / MZ2;
  const double VD   = - 0.5 + 2 * s2tw / 3;
  const double VU   = + 0.5 - 4 * s2tw / 3;
  const double AD   = - 0.5;
  const double AU   = + 0.5;
  const double Ve   = - 0.5 + 2 * s2tw;
  const double Ae   = - 0.5;
  
  auto &x = *GetBinValues(td, "x");
  auto &Q2 = *GetBinValues(td, "Q2");
  //auto &y = *GetBinValues(td, "y");
  std::vector<double> StrFct(12);
  for (size_t i = 0; i <= x.size(); i++)
  {
    const double Q = std::sqrt(Q2[i]);
    const double PZ  = Q / ( Q + MZ2 ) / ( 4 * s2tw * ( 1 - s2tw ) );
    const double PZ2 = PZ * PZ;
    hoppetStrFct(x[i],Q, xmuR * Q, xmuF * Q, &StrFct[0]);

//  const double F1NCh = StrFct[iF1EM] + StrFct[iF1Z] * (Ve * Ve + Ae * Ae) * PZ2 - StrFct[iF1gZ] * Ve * PZ;
    const double F2NCh = StrFct[iF2EM] + StrFct[iF2Z] * (Ve * Ve + Ae * Ae) * PZ2 - StrFct[iF2gZ] * Ve * PZ;
  // const double F3NCh = 2 * StrFct[iF3Z] * Ae * Ve * PZ2 - StrFct[iF3gZ] * Ae * PZ;
  // const double FLNCh = F2NCh - 2 * x * F1NCh;
  
    valExternal[i] = F2NCh;
  }
}


void ReactionHOPPET_DISNC::atIteration()
{
  // Make sure to call the parent class initialization:
  super::atIteration();
}
