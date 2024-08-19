/*
   @file Reaction_DISNC_Hoppet.cc
   @date 2017-04-08
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-04-08
*/
#include "ReactionTheory.h"
#include "ReactionHOPPET_DISNC.h"
#include <iostream>
#include <cstdio>
#include "ReactionBaseDISNC.h"
#include "hoppet_v1.h" // Include the HOPPET header
#include "xfitter_pars.h"
#include <BaseEvolution.h>
#include <hf_errlog.h>

using namespace hoppetv1;
using namespace std;

// The class factories
extern "C" ReactionHOPPET_DISNC *create()
{
    return new ReactionHOPPET_DISNC();
}

void lha_unpolarized_dummy_pdf(const double &x, const double &Q, double *pdf)
{
    double uv, dv;
    double ubar, dbar;
    double N_g = 1.7, N_ls = 0.387975;
    double N_uv = 5.107200, N_dv = 3.064320;
    double N_db = N_ls / 2;

    uv = N_uv * pow(x, 0.8) * pow((1 - x), 3);
    dv = N_dv * pow(x, 0.8) * pow((1 - x), 4);
    dbar = N_db * pow(x, -0.1) * pow(1 - x, 6);
    ubar = dbar * (1 - x);

    pdf[0 + 6] = N_g * pow(x, -0.1) * pow(1 - x, 5);
    pdf[-3 + 6] = 0.2 * (dbar + ubar);
    pdf[3 + 6] = pdf[-3 + 6];
    pdf[2 + 6] = uv + ubar;
    pdf[-2 + 6] = ubar;
    pdf[1 + 6] = dv + dbar;
    pdf[-1 + 6] = dbar;

    pdf[4 + 6] = 0;
    pdf[5 + 6] = 0;
    pdf[6 + 6] = 0;
    pdf[-4 + 6] = 0;
    pdf[-5 + 6] = 0;
    pdf[-6 + 6] = 0;
}

// Main function to compute results at an iteration
void ReactionHOPPET_DISNC::F2(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal)
{
    const double mc = *XFITTER_PARS::getParamD("mch");
    const double mb = *XFITTER_PARS::getParamD("mbt");
    const double mt = *XFITTER_PARS::getParamD("mtp");
    cout<<"61"<<endl;
    
    hoppetSetPoleMassVFN(mc, mb, mt);

    const bool param_coefs = true;
    const double xmuR = 1.;
    const double xmuF = 1.;

    double Qmax = 13000.0;
    double Qmin = 1.0;
    int order = -6;
    double ymax = 16.0;
    double dy = 0.05;
    double dlnlnQ = dy / 4.0;
    int nloop = 3;
    double minQval = min(xmuF * Qmin, Qmin);
    double maxQval = max(xmuF * Qmax, Qmax);

    // Initialize HOPPET
   // hoppetStartExtended(ymax, dy, minQval, maxQval, dlnlnQ, nloop, order, factscheme_MSbar);
    
    int nflav = -5;
    int order_max = 4;
    int sc_choice = scale_choice_Q;
    double zmass = 91.1876;
    double wmass = 80.377;
    
   // hoppetStartStrFctExtended(order_max, nflav, scale_choice_Q, zmass, param_coefs, wmass, zmass);

    double asQ = 0.35;
    double Q0 = sqrt(2.0);
    double muR_Q = 1.0;

    // Evolve the PDF
  //  hoppetEvolve(asQ, Q0, nloop, muR_Q, lha_unpolarized_dummy_pdf, Q0);

    // Initialize structure functions
  //  hoppetInitStrFct(order_max, param_coefs, xmuR, xmuF);

    // Obtain parameters
    const auto _convfac = *XFITTER_PARS::getParamD("convFac");
    const auto _alphaem = *XFITTER_PARS::getParamD("alphaem");
    const auto MZ = *XFITTER_PARS::getParamD("Mz");
    const auto MW = *XFITTER_PARS::getParamD("Mw");

    const double MW2 = MW * MW;
    const double MZ2 = MZ * MZ;

    // Construct structure functions
    const double s2tw = 1 - MW2 / MZ2;
    const double VD = -0.5 + 2 * s2tw / 3;
    const double VU = +0.5 - 4 * s2tw / 3;
    const double AD = -0.5;
    const double AU = +0.5;
    const double Ve = -0.5 + 2 * s2tw;
    const double Ae = -0.5;

    auto &xvals = *GetBinValues(td, "x");
    auto &Qvals = *GetBinValues(td, "Q2");
    double pdf[13];
    std::vector<double> StrFct(14); // Adjusted to match example
    cout<<Qvals.size()<<endl;

    printf("  x        Q       F2NCh       F3NCh       FLNCh\n");
    printf("------------------------------------------------\n");

    for (size_t i = 0; i < xvals.size(); i++)
    {
		
        const double x = xvals[i];
         cout<<"x="<<x<<endl;
         cout<<Qvals.size()<<endl;
        const double Q = (sqrt(Qvals[i]));
        cout<<Qvals.size()<<endl;
  
        // Compute propagator factors
        const double Q2 = Q * Q;
        const double PZ = Q2 / (Q2 + MZ2) / (4.0 * s2tw * (1.0 - s2tw));
        const double PZ2 = PZ * PZ;
         
        // cout<<"141"<<endl;
        // Evaluate PDFs at (x, Q)
       // hoppetEval(x, Q, pdf);
         cout<<"144"<<endl;
        // Evaluate structure functions at (x, Q)
        hoppetStrFct(x, Q, Q, Q, &StrFct[0]);
         cout<<"147"<<endl;
        // Extract individual structure functions
        const double F1EM = StrFct[iF1EM];
        const double F2EM = StrFct[iF2EM];
        const double F3EM = 0.0;           // Electromagnetic interactions are parity conserving
        const double F1Z = StrFct[iF1Z];
        const double F2Z = StrFct[iF2Z];
        const double F3Z = StrFct[iF3Z];
        const double F1gZ = StrFct[iF1gZ];
        const double F2gZ = StrFct[iF2gZ];
        const double F3gZ = StrFct[iF3gZ];

        // Compute combined structure functions
        const double F1NCh = F1EM + F1Z * (Ve * Ve + Ae * Ae) * PZ2 - F1gZ * Ve * PZ;
        const double F2NCh = F2EM + F2Z * (Ve * Ve + Ae * Ae) * PZ2 - F2gZ * Ve * PZ;
        const double F3NCh = 2.0 * F3Z * Ae * Ve * PZ2 - F3gZ * Ae * PZ;
        const double FLNCh = F2NCh - 2.0 * x * F1NCh;

        // Store F2NCh as the output
        valExternal[i] = F2NCh + F3NCh + FLNCh;

        // Print the results		
        printf("%7.5f  %7.3f  %10.6f  %10.6f  %10.6f  %10.6f\n",
               x, Q, F1NCh, F2NCh, F3NCh, FLNCh);
    }
}

void ReactionHOPPET_DISNC::atIteration()
{
    // Make sure to call the parent class initialization:
    super::atIteration();
}
