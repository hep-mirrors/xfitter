/*
   @file Reaction_DISCC_HOPPET.cc
   @date 2017-04-08
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-04-08
*/
#include "ReactionTheory.h"
#include "ReactionHOPPET_DISCC.h"
#include <iostream>
#include <cstdio>
#include "ReactionBaseDISCC.h"
#include "hoppet_v1.h" // Include the HOPPET header
#include "xfitter_pars.h"
#include <BaseEvolution.h>
#include <hf_errlog.h>
#include "xfitter_cpp_base.h"

using namespace hoppetv1;
using namespace std;

// The class factories
extern "C" ReactionHOPPET_DISCC *create()
{
    return new ReactionHOPPET_DISCC();
}

// Initialize at the start of the computation
void ReactionHOPPET_DISCC::atStart()
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

    const bool exact_nfthreshold = true;
    const bool exact_splitting = false;
    hoppetSetExactDGLAP(exact_nfthreshold,exact_splitting);

    // temporary, needed for alphaS evolution
    _alphas = XFITTER_PARS::getParamD("alphas");
    _Q0 = *XFITTER_PARS::getParamD("Q0");
}

void ReactionHOPPET_DISCC::initTerm(TermData *td)
{
    Super::initTerm(td);
    BaseDISCC::ReactionData *rd = new BaseDISCC::ReactionData();
    unsigned termID = td->id;
    if(rd->_dataFlav != BaseDISCC::dataFlav::incl) {
        hf_errlog(20240905, "F: HOPPET supports only inclusive DIS structure functions");
    }
    auto nBins = td->getNbins();
    if(rd->_integrated) {
        nBins = rd->_integrated->getBinValuesQ2()->size();
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

void ReactionHOPPET_DISCC::atIteration() {
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
        _order_HOPPET_Evolution = _order;
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

valarray<double> ReactionHOPPET_DISCC::F2(TermData *td)
{
    //for (size_t i = 0; i < GetBinValues(td, "x")->size(); i++) {valExternal[i] = 0.;}
    calcF2FLF3(td->id);
    return _f2[td->id];
}
valarray<double> ReactionHOPPET_DISCC::FL(TermData *td) {
    //for (size_t i = 0; i < GetBinValues(td, "x")->size(); i++) {valExternal[i] = 0.;}
    calcF2FLF3(td->id);
    return _fl[td->id];
}
valarray<double> ReactionHOPPET_DISCC::xF3(TermData *td) {
    //for (size_t i = 0; i < GetBinValues(td, "x")->size(); i++) {valExternal[i] = 0.;}
    calcF2FLF3(td->id);
    return _f3[td->id];
}

void ReactionHOPPET_DISCC::calcF2FLF3(unsigned dataSetID) {
    // skip if already calculated, e.g. for calls from FL(), xF3() after F2()
    if ((_f2[dataSetID][0] > -99.)) {
        return;
    }

    // setup HOPPET evolution - try to copy it from our evolution
    auto td = _tdDS[dataSetID];
    td->actualizeWrappers();
    BaseDISCC::ReactionData *rd = new BaseDISCC::ReactionData();
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

    const int charge = rd->_charge;

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
        switch (rd->_dataFlav)
        {
        case BaseDISCC::dataFlav::incl:
            if (charge == 1) {
                _f2[dataSetID][i] = StrFct[iF2Wp];
                _fl[dataSetID][i] = StrFct[iF2Wp] - 2 * x * StrFct[iF1Wp];
                _f3[dataSetID][i] = StrFct[iF3Wp] * x;
            }
            else {
                _f2[dataSetID][i] = StrFct[iF2Wm];
                _fl[dataSetID][i] = StrFct[iF2Wm] - 2 * x * StrFct[iF1Wm];
                _f3[dataSetID][i] = StrFct[iF3Wm] * x;
            }
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
