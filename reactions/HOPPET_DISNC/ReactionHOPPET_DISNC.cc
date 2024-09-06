/*
   @file Reaction_DISNC_HOPPET.cc
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

// Initialize at the start of the computation
void ReactionHOPPET_DISNC::atStart()
{
    // do not call parent atStart(): it initialises QCDNUM
    // Super::atStart();

    _order = OrderMap(XFITTER_PARS::getParamS("Order"));
    _convfac = XFITTER_PARS::getParamD("convFac");
    _alphaem = XFITTER_PARS::getParamD("alphaem");
    _Mz = XFITTER_PARS::getParamD("Mz");
    _Mw = XFITTER_PARS::getParamD("Mw");
    _sin2thetaW = XFITTER_PARS::getParamD("sin2thW");

    const double mc = *XFITTER_PARS::getParamD("mch");
    const double mb = *XFITTER_PARS::getParamD("mbt");
    const double mt = *XFITTER_PARS::getParamD("mtp");
    hoppetSetPoleMassVFN(mc, mb, mt);
    //hoppetSetVFN(mc, mb, mt);

    // temporary, needed for alphaS evolution
    _alphas = XFITTER_PARS::getParamD("alphas");
    _Q0 = *XFITTER_PARS::getParamD("Q0");
}

void ReactionHOPPET_DISNC::initTerm(TermData *td)
{
    Super::initTerm(td);
    unsigned termID = td->id;
    if(GetDataFlav(termID) != dataFlav::incl) {
        hf_errlog(20240905, "F: HOPPET supports only inclusive DIS structure functions");
    }
    auto nBins = td->getNbins();
    if(_integrated.find(termID) != _integrated.end()) {
        nBins = _integrated[termID]->getBinValuesQ2()->size();
    }
    _f2[termID].resize(nBins);
    _fl[termID].resize(nBins);
    _f3[termID].resize(nBins);

    // read reaction specific parameters
    _dy = *td->getParamD("dy");
    _xmuR = *td->getParamD("xmuR");
    _xmuF = *td->getParamD("xmuF");
    _muR_Q = *td->getParamD("muR_Q");
    _param_coefs = td->getParamI("param_coefs");
}

void ReactionHOPPET_DISNC::atIteration() {
    // Make sure to call the parent class initialization:
    super::atIteration();

    // Reset internal arrays
    for (auto ds : _dsIDs)
    {
        (_f2[ds])[0] = -100.;
        (_fl[ds])[0] = -100.;
        (_f3[ds])[0] = -100.;
    }

    // the code below needs to be done only once, however,
    // it requires reaction specific parameters which can 
    // be read only by initTerm() and are not available atStart() 
    // -> therefore call it only once here (this works correctly 
    // if these parameters are not changed at iteration)
    static int init = 0;
    if (init == 0) {
        // temporary: allow different orders in evolution and DIS SFs
        if (XFITTER_PARS::gParametersS.find("Order_HOPPET_Evolution") != XFITTER_PARS::gParametersS.end()) {
            _order_HOPPET_Evolution = OrderMap(XFITTER_PARS::getParamS("Order_HOPPET_Evolution"));
        }

        hoppetStart(_dy, _order_HOPPET_Evolution);

        // Extended HOPPET Initialization (from example): seems to be not needed
        //double Qmax = 13000.0;
        //double Qmin = 1.0;
        //int order = -6;
        //double ymax = 16.0;
        //double dlnlnQ = dy / 4.0;
        //int nloop = 3;
        //double xmuF   = 1.0;
        //double minQval = min(xmuF * Qmin, Qmin);
        //double maxQval = max(xmuF * Qmax, Qmax);
        //hoppetStartExtended(ymax, dy, minQval, maxQval, dlnlnQ, nloop, order, factscheme_MSbar);
        
        int nflav = -1 * XFITTER_PARS::getParamI("NFlavour"); // negative nflav to use a variable-flavour number scheme
        hoppetStartStrFctExtended(_order, nflav, scale_choice_Q, *_Mz, _param_coefs, *_Mw, *_Mz);

        init = 1;
    }
}

void ReactionHOPPET_DISNC::F2(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal)
{
    //for (size_t i = 0; i < GetBinValues(td, "x")->size(); i++) {valExternal[i] = 0.;}
    calcF2FLF3(td->id);
    valExternal = _f2[td->id];
}
void ReactionHOPPET_DISNC::FL(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal) {
    //for (size_t i = 0; i < GetBinValues(td, "x")->size(); i++) {valExternal[i] = 0.;}
    calcF2FLF3(td->id);
    valExternal = _fl[td->id];
}
void ReactionHOPPET_DISNC::xF3(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal) {
    //for (size_t i = 0; i < GetBinValues(td, "x")->size(); i++) {valExternal[i] = 0.;}
    calcF2FLF3(td->id);
    valExternal = _f3[td->id];
}

void ReactionHOPPET_DISNC::calcF2FLF3(unsigned dataSetID) {
    // skip if already calculated, e.g. for calls from FL(), xF3() after F2()
    if ((_f2[dataSetID][0] > -99.)) {
        return;
    }

    // setup HOPPET evolution - try to copy it from our evolution
    auto td = _tdDS[dataSetID];
    td->actualizeWrappers();
    // TODO: it seems that still one needs to call hoppetEvolve() in order to get alphaS evolution
    // how to get alphaS assigned via hoppetAssign()?
    const double Q_for_alphaS = *_Mz;
    hoppetEvolve( *_alphas, Q_for_alphaS, _order_HOPPET_Evolution, _muR_Q, pdf_xfxq_wrapper1_, _Q0);
    //double f[13];
    //hoppetEval(0.001,10.,f);
    //printf("f[6] = %f  as = %f\n", f[6], hoppetAlphaS(10.));
    hoppetAssign(pdf_xfxq_wrapper1_);
    //printf("HOPPET assigned\n");
    //hoppetEval(0.001,10.,f);
    //printf("f[6] = %f  as = %f\n", f[6], hoppetAlphaS(10.));

    // Initialize structure functions
    hoppetInitStrFct(_order, _param_coefs, _xmuR, _xmuF);

    const int charge = GetCharge(dataSetID);
    const double pol = GetPolarisation(dataSetID);

    const double _ve = -0.5 + 2. * (*_sin2thetaW); // !
    const double _ae = -0.5;                    // !
    const double _au = 0.5;
    const double _ad = -0.5;
    const double _vu = _au - (4. / 3.) * (*_sin2thetaW);
    const double _vd = _ad + (2. / 3.) * (*_sin2thetaW);
    const double cos2thetaW = 1 - *_sin2thetaW;
    const double MZ2 = *_Mz * *_Mz;

    // Initialize structure functions
    hoppetInitStrFct(_order, _param_coefs, _xmuR, _xmuF);

    auto &xvals = *GetBinValues(td, "x");
    auto &Q2vals = *GetBinValues(td, "Q2");
    double StrFct[14];

    for (size_t i = 0; i < xvals.size(); i++)
    {
		
        const double x = xvals[i];
        const double Q2 = Q2vals[i];
        const double Q = sqrt(Q2);
        if(Q < 1.) {
            _f2[dataSetID][i] = _fl[dataSetID][i] = _f3[dataSetID][i] = 9e99;
            continue;
        }

        hoppetStrFct(x, Q, _xmuR * Q, _xmuF * Q, &StrFct[0]);
        //printf("StrFct = ");
        //for(size_t ii = 0; ii < 14; ii++) printf("  %f", StrFct[ii]);
        //printf("\n");
        double k = 1. / (4 * *_sin2thetaW * cos2thetaW) * Q2 / (Q2 + MZ2);
        switch (GetDataFlav(td->id))
        {
        case dataFlav::incl:
            _f2[dataSetID][i] = StrFct[iF2EM] - (_ve + charge * pol * _ae) * k * StrFct[iF2gZ] + (_ae * _ae + _ve * _ve + 2 * charge * pol * _ae * _ve) * k * k * StrFct[iF2Z];
            _fl[dataSetID][i] = _f2[dataSetID][i] - 2 * x * (StrFct[iF1EM] - (_ve + charge * pol * _ae) * k * StrFct[iF1gZ] + (_ae * _ae + _ve * _ve + 2 * charge * pol * _ae * _ve) * k * k * StrFct[iF1Z]);
            _f3[dataSetID][i] = (_ae * charge + pol * _ve) * k * x * StrFct[iF3gZ] + (-2 * _ae * _ve * charge - pol * (_ve * _ve + _ae * _ae)) * k * k * x * StrFct[iF3Z];
            break;
        // should not be here
        //case dataFlav::c:
        //    _f2[dataSetID][i] = _fl[dataSetID][i] = _f3[dataSetID][i] = 0.;
        //    break;
        //case dataFlav::b:
        //    _f2[dataSetID][i] = _fl[dataSetID][i] = _f3[dataSetID][i] = 0.;
        //    break;
        }
    }
}
