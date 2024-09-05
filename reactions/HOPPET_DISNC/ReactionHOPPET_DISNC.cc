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
#include "xfitter_cpp_base.h"

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

void example() {
  double mc = 1.414213563;   
  double mb = 4.5;
  double mt = 175.0;
  
  hoppetSetPoleMassVFN(mc,mb,mt);

  double Qmax   = 13000.0;
  double Qmin   = 1.0;
  int order      = -6; 
  double ymax    = 16.0;
  double dy      = 0.05;
  double dlnlnQ  = dy/4.0;
  int    nloop   = 3;
  double xmuR   = 1.0;
  double xmuF   = 1.0;
  double minQval = min(xmuF*Qmin, Qmin);
  double maxQval = max(xmuF*Qmax, Qmax);

  hoppetStartExtended(ymax,dy,minQval,maxQval,dlnlnQ,nloop,order,factscheme_MSbar);
  
  int nflav = -5;
  int order_max = 4;
  int sc_choice = scale_choice_Q;
  double zmass = 91.1876;
  double wmass = 80.377;
  bool param_coefs = true;

  hoppetStartStrFctExtended(order_max, nflav,sc_choice,zmass,param_coefs,wmass,zmass);
    
  //double asQ      = 0.35;
  //double Q0       = sqrt(2.0);
  double asQ      = 0.118;
  double Q0       = zmass;
  double muR_Q    = 1.0;

  double _convfac = *XFITTER_PARS::getParamD("convFac");
  double _alphaem = *XFITTER_PARS::getParamD("alphaem");
  double _Mz = *XFITTER_PARS::getParamD("Mz");
  double _Mw = *XFITTER_PARS::getParamD("Mw");
  double _sin2thetaW = *XFITTER_PARS::getParamD("sin2thW");

  double _ve = -0.5 + 2. * _sin2thetaW; // !
  double _ae = -0.5;                    // !
  double _au = 0.5;
  double _ad = -0.5;
  double _vu = _au - (4. / 3.) * _sin2thetaW;
  double _vd = _ad + (2. / 3.) * _sin2thetaW;

  hoppetEvolve(asQ, Q0, nloop, muR_Q, lha_unpolarized_dummy_pdf, Q0);
  hoppetAssign(pdf_xfxq_wrapper_);

  hoppetInitStrFct(order_max,param_coefs, xmuR,xmuF);
    
  // output the results
  double pdf[13];
  //double xvals[9]={1e-5,1e-4,1e-3,1e-2,0.1,0.3,0.5,0.7,0.9};
  double xvals[2]={0.4, 0.65};
  //double Q = 100;
  double Q = 173.205;
  double StrFct[14];
  printf("                                Evaluating PDFs and structure functions at Q = %8.3f GeV\n",Q);
  printf("    x      u-ubar      d-dbar    2(ubr+dbr)    c+cbar       gluon       F1γ         F2γ         F1Z         F2Z         F3Z -->> val\n");
  //for (int ix = 0; ix < 9; ix++) {
  for (int ix = 0; ix < 2; ix++) {
    hoppetEval(xvals[ix], Q, pdf);
    hoppetStrFct(xvals[ix],Q,Q,Q,StrFct);
    double f2g = StrFct[iF2EM];
    double f2gZ = StrFct[iF2gZ ];
    double f2Z = StrFct[iF2Z ];
    int charge = +1;
    double pol = 0.;
    double cos2thetaW = 1 - _sin2thetaW;
    double k = 1. / (4 * _sin2thetaW * cos2thetaW) * (Q*Q) / ((Q*Q) + _Mz * _Mz);
    double val = f2g - (_ve + charge * pol * _ae) * k * f2gZ + (_ae * _ae + _ve * _ve + 2 * charge * pol * _ae * _ve) * k * k * f2Z;
    printf("%7.1E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E -->> %11.4E\n",xvals[ix],
           pdf[6+2]-pdf[6-2], 
           pdf[6+1]-pdf[6-1], 
           2*(pdf[6-1]+pdf[6-2]),
           (pdf[6-4]+pdf[6+4]),
           pdf[6+0],
	   StrFct[iF1EM],
	   StrFct[iF2EM],
	   StrFct[iF1Z],
	   StrFct[iF2Z ],
	   StrFct[iF3Z ],
       val);
  }
}


// Initialize at the start of the computation
// (copied from FFABM)
void ReactionHOPPET_DISNC::atStart()
{
  // do not call parent atStart(): it initialises QCDNUM
  // Super::atStart();
}

void ReactionHOPPET_DISNC::initTerm(TermData *td)
{
  Super::initTerm(td);

  unsigned termID = td->id;
  auto nBins = td->getNbins();
  if(_integrated.find(termID) != _integrated.end()) {
    nBins = _integrated[termID]->getBinValuesQ2()->size();
  }
  _f2[termID].resize(nBins);
  _fl[termID].resize(nBins);
  _f3[termID].resize(nBins);
}

void ReactionHOPPET_DISNC::atIteration()
{
    // Make sure to call the parent class initialization:
    super::atIteration();

    // Reset internal arrays
    for (auto ds : _dsIDs)
    {
        (_f2[ds])[0] = -100.;
        (_fl[ds])[0] = -100.;
        (_f3[ds])[0] = -100.;
    }
}

void ReactionHOPPET_DISNC::F2(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal)
{
    //for (size_t i = 0; i < GetBinValues(td, "x")->size(); i++) {valExternal[i] = 0.;}
    calcF2FLF3(td->id);
    valExternal = _f2[td->id];
    printf("F2 valExternal = %f %f\n", valExternal[0], valExternal[1]);
}
void ReactionHOPPET_DISNC::FL(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal) {
    //for (size_t i = 0; i < GetBinValues(td, "x")->size(); i++) {valExternal[i] = 0.;}
    calcF2FLF3(td->id);
    valExternal = _fl[td->id];
    printf("FL valExternal = %f %f\n", valExternal[0], valExternal[1]);
}
void ReactionHOPPET_DISNC::xF3(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal) {
    //for (size_t i = 0; i < GetBinValues(td, "x")->size(); i++) {valExternal[i] = 0.;}
    calcF2FLF3(td->id);
    valExternal = _f3[td->id];
    printf("F3 valExternal = %f %f\n", valExternal[0], valExternal[1]);
}

void ReactionHOPPET_DISNC::calcF2FLF3(unsigned dataSetID) {
    // skip if already calculated, e.g. for calls from FL(), xF3() after F2()
    if ((_f2[dataSetID][0] > -99.)) {
        return;
    }

    // otherwise do actual calculations
    auto td = _tdDS[dataSetID];
    td->actualizeWrappers();

    //example();
    //return;

    const double _convfac = *XFITTER_PARS::getParamD("convFac");
    const double _alphaem = *XFITTER_PARS::getParamD("alphaem");
    const double _Mz = *XFITTER_PARS::getParamD("Mz");
    const double _Mw = *XFITTER_PARS::getParamD("Mw");
    const double _sin2thetaW = *XFITTER_PARS::getParamD("sin2thW");
    const int charge = GetCharge(td->id);
    const double pol = GetPolarisation(td->id);

    const double _ve = -0.5 + 2. * _sin2thetaW; // !
    const double _ae = -0.5;                    // !
    const double _au = 0.5;
    const double _ad = -0.5;
    const double _vu = _au - (4. / 3.) * _sin2thetaW;
    const double _vd = _ad + (2. / 3.) * _sin2thetaW;
    const double cos2thetaW = 1 - _sin2thetaW;

    //const int PtOrder = OrderMap(XFITTER_PARS::getParamS("Order")) - 1;//here was -1
    const int PtOrder = OrderMap(XFITTER_PARS::getParamS("Order"));//here was -1
    const double mc = *XFITTER_PARS::getParamD("mch");
    const double mb = *XFITTER_PARS::getParamD("mbt");
    const double mt = *XFITTER_PARS::getParamD("mtp");
    hoppetSetPoleMassVFN(mc, mb, mt);

    //double dy = 0.05;
    double dy = 0.1;
    hoppetStart(dy, PtOrder);
    //hoppetSetVFN(mc, mb, mt);

    double Qmax = 13000.0;
    double Qmin = 1.0;
    int order = -6;
    double ymax = 16.0;
    //double ymax = 12.0; // this works
    double dlnlnQ = dy / 4.0;
    int nloop = 3;
    double xmuF   = 1.0;
    double minQval = min(xmuF * Qmin, Qmin);
    double maxQval = max(xmuF * Qmax, Qmax);
    // Initialize HOPPET
    //hoppetStartExtended(ymax, dy, minQval, maxQval, dlnlnQ, nloop, order, factscheme_MSbar);
    
    int nflav = -5; // negative nflav to use a variable-flavour number scheme
    //int sc_choice = scale_choice_Q;
    //const double zmass = *XFITTER_PARS::getParamD("Mz");
    //const double wmass = *XFITTER_PARS::getParamD("Ww");
    const bool param_coefs = true;
    const double xmuR = 1.;
    //const double xmuF = 1.;    
    hoppetStartStrFctExtended(PtOrder, nflav, scale_choice_Q, _Mz, param_coefs, _Mw, _Mz);

    //double asQ = 0.35;
    //double Q0 = sqrt(2.0);
    //double muR_Q = 1.0;

    // Evolve the PDF
    const double* Alphas_ref = XFITTER_PARS::getParamD("alphas");
    const double* Q0         = XFITTER_PARS::getParamD("Q0");
    const double* Q_ref      = XFITTER_PARS::getParamD("Mz");
    //hoppetEvolve( *Alphas_ref, *Q_ref, PtOrder, 1.0, heralhc_init, *Q0);
    hoppetEvolve( *Alphas_ref, *Q_ref, PtOrder, 1.0, lha_unpolarized_dummy_pdf, *Q0);
    //hoppetEvolve(asQ, Q0, nloop, muR_Q, lha_unpolarized_dummy_pdf, Q0);
    double f[13];
    hoppetEval(0.001,10.,f);
    printf("f[6] = %f  as = %f\n", f[6], hoppetAlphaS(10.));
    hoppetAssign(pdf_xfxq_wrapper_);
    printf("HOPPET assigned\n");
    hoppetEval(0.001,10.,f);
    printf("f[6] = %f  as = %f\n", f[6], hoppetAlphaS(10.));

    // Initialize structure functions
    hoppetInitStrFct(PtOrder, param_coefs, xmuR, xmuF);

    auto &xvals = *GetBinValues(td, "x");
    auto &Q2vals = *GetBinValues(td, "Q2");
    double pdf[13];
    //std::vector<double> StrFct(14); // Adjusted to match example
    double StrFct[14];
    //cout<<Qvals.size()<<endl;

    printf("  x        Q       F2NCh       F3NCh       FLNCh\n");
    printf("------------------------------------------------\n");

    for (size_t i = 0; i < xvals.size(); i++)
    {
		
        const double x = xvals[i];
         //cout<<Q2vals.size()<<endl;
        const double Q = sqrt(Q2vals[i]);
         //cout<<"x="<<x<<" Q="<<Q<<endl;
        if(Q<1.) {
            _f2[dataSetID][i] = _f2[dataSetID][i] = _f2[dataSetID][i] = 9e99;
            continue;
        }
        //if(Q<10.) continue;
        //cout<<Q2vals.size()<<endl;
  
        // Compute propagator factors
        //const double Q2 = Q * Q;
        //const double PZ = Q2 / (Q2 + MZ2) / (4.0 * s2tw * (1.0 - s2tw));
        //const double PZ2 = PZ * PZ;
         
        // cout<<"141"<<endl;
        // Evaluate PDFs at (x, Q)
     //  hoppetEval(x, Q, pdf);
       //  cout<<"144"<<endl;
        // Evaluate structure functions at (x, Q)
        hoppetStrFct(x, Q, Q, Q, &StrFct[0]);
        printf("StrFct = ");
        for(size_t ii = 0; ii < 14; ii++) {
            printf("  %f", StrFct[ii]);
        }
        printf("\n");
        // cout<<"147"<<endl;
        // Extract individual structure functions
        /*const double F1EM = StrFct[iF1EM];
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
        const double FLNCh = F2NCh - 2.0 * x * F1NCh;*/

        // Store F2NCh as the output
        double k = 1. / (4 * _sin2thetaW * cos2thetaW) * (Q*Q) / ((Q*Q) + _Mz * _Mz);
        switch (GetDataFlav(td->id))
        {
        case dataFlav::incl:
            // f2g - (_ve + charge * pol * _ae) * k * f2gZ + (_ae * _ae + _ve * _ve + 2 * charge * pol * _ae * _ve) * k * k * f2Z;
            _f2[dataSetID][i] = StrFct[iF2EM] - (_ve + charge * pol * _ae) * k * StrFct[iF2gZ] + (_ae * _ae + _ve * _ve + 2 * charge * pol * _ae * _ve) * k * k * StrFct[iF2Z];
            _fl[dataSetID][i] = _f2[dataSetID][i] - 2 * x * (StrFct[iF1EM] - (_ve + charge * pol * _ae) * k * StrFct[iF1gZ] + (_ae * _ae + _ve * _ve + 2 * charge * pol * _ae * _ve) * k * k * StrFct[iF1Z]);
            _f3[dataSetID][i] = (_ae * charge + pol * _ve) * k * x * StrFct[iF3gZ] + (-2 * _ae * _ve * charge - pol * (_ve * _ve + _ae * _ae)) * k * k * x * StrFct[iF3Z];
            break;
        case dataFlav::c:
            _f2[dataSetID][i] = _fl[dataSetID][i] = _f3[dataSetID][i] = 0.;
            break;
        case dataFlav::b:
            _f2[dataSetID][i] = _fl[dataSetID][i] = _f3[dataSetID][i] = 0.;
            break;
        }
        //_f2[dataSetID][i] = F2NCh;
        //_fl[dataSetID][i] = FLNCh;
        //_f3[dataSetID][i] = F3NCh*x;
        //_f3[dataSetID][i] = 0.;
        //_fl[dataSetID][i] = 0.;
        //valExternal[i] = F2NCh + F3NCh + FLNCh;

        // Print the results		
        //printf("%7.5f  %7.3f  %10.6f  %10.6f  %10.6f  %10.6f  -->>  %f\n", x, Q, F1NCh, F2NCh, F3NCh, FLNCh, valExternal[i]);
        //printf("%7.5f  %7.3f  %10.6f  %10.6f  %10.6f  %10.6f  -->>  %f\n", x, Q, F1NCh, F2NCh, F3NCh, FLNCh, F2NCh + F3NCh + FLNCh);
    }
}
