/*
  @file ReactionPineAPPL.cc
  @date 2017-03-28
  @author  AddReaction.py
  Created by  AddReaction.py on 2017-03-28
*/

#include "ReactionPineAPPL.h"
#include "pineappl_capi.h"
#include "xfitter_pars.h"
#include "xfitter_steer.h"
#include "xfitter_cpp_base.h"
#include <memory>
#include "BaseEvolution.h"
#include "TermData.h"

using namespace std;
using namespace xfitter;

struct DatasetData {
    vector<pineappl_grid*> grids;
    int order;
    int32_t pdg_id;
    const double *muR, *muF; // !> renormalisation and factorisation scales
    bool flagNorm;  // !> if true, multiply by bin width
};

const double ONE=1;

extern "C" ReactionPineAPPL* create() {
    return new ReactionPineAPPL();
}
//MESSING STARTS HERE
void ReactionPineAPPL::initTerm(TermData*td) {
    DatasetData* data = new DatasetData;
    td->reactionData = (void*)data;
    //Load grids
    string GridName = td->getParamS("GridName");
    try {
        istringstream ss(GridName);
        string token;
        while (getline(ss, token, ',')) {
             pineappl_grid* g = pineappl_grid_read(token.c_str());             
             data->grids.push_back(g);
        }
    } catch ( const exception& e ) {
        hf_errlog(22121301, "E: Unhandled exception while trying to read PineAPPL grid(s) "+GridName);
        throw e;
    }

    if (td->hasParam("PDG")) data->pdg_id = td->getParamI("PDG");
    else {
        data->pdg_id = 2212;
        hf_errlog(22121201, "W: PineAPPL initTerm assuming protons");
    }
    // Get Order
    if (td->hasParam("Order")) {  //FIXME non-functional; using all orders in grid, cf below
        string ordS = td->getParamS("Order");
        if (ordS.find("LO")!=string::npos) data->order = ordS.erase(ordS.find("LO")).size();
        else data->order = -1; //Use all orders available in grid
    }
    // Get MuR and MuF, or use 1.0 as default
    if (td->hasParam("muR")) data->muR = td->getParamD("muR");
    else data->muR = &ONE;
    if (td->hasParam("muF")) data->muF = td->getParamD("muF");
    else data->muF = &ONE;
    // Get if should normalize by dividing by bin width (no by default)
    data->flagNorm = false;
    if (td->hasParam("norm")) {
        int norm=td->getParamI("norm");
        if (norm==1) data->flagNorm=true;
        else if (norm!=0) hf_errlog(22121202, "F: unrecognised norm = " + norm);
    }
    size_t Ngrids = data->grids.size();
} //initTerm

void ReactionPineAPPL::freeTerm(TermData*td) {
    DatasetData* data = (DatasetData*)td->reactionData;
    size_t Ngrids=data->grids.size();
    for (size_t i=0; i<Ngrids; ++i) pineappl_grid_delete(data->grids[i]);
    delete data;
}

//Pineappl assumes that the PDF and alphaS wrapper function pointers should be given a pointer
//"state" of e.g. an LHAPDF object if the values were read from there. Here we just call the 
//existing wrappers, but rearrange the parameter list and return value to suit the demands of
//the pineappl convolution function.
double xfx(int32_t id, double x, double q2, void *state) {
    double pdfs[13];
    pdf_xfxq_wrapper_(x, sqrt(q2), pdfs);
    return pdfs[id+6];
}
double alphas(double q2, void *state) {
    return alphas_wrapper_(sqrt(q2));
}

void ReactionPineAPPL::compute(TermData*td,valarray<double>&val,map<string,valarray<double> >&err) {
    const DatasetData& data = *(DatasetData*)td->reactionData;
    const int order = data.order;
    const double muR = *data.muR;
    const double muF = *data.muF;
    unsigned int pos = 0;

    //calculate output array size
    size_t np=0;
    for (pineappl_grid* grid : data.grids) if (grid) np += pineappl_grid_bin_count(grid);
    val.resize(np);

    for (pineappl_grid* grid : data.grids) {
        vector<double> gridVals;
        if (grid) {//real, non-dummy grid
            gridVals.resize(pineappl_grid_bin_count(grid));
            td->actualizeWrappers();
            //See function specification in deps/pineappl/include/pineappl_capi/pineappl_capi.h
            pineappl_grid_convolute_with_one(grid, data.pdg_id, 
                                             &xfx, &alphas, 
                                             nullptr,//state: same as provided to wrappers. Redundant for xFitter
                                             nullptr,//order_mask FIXME simply uses all in grid if null
                                             nullptr,//lumi mask  FIXME simply uses all in grid if null
                                             muR, muF, gridVals.data());
          //scale by bin width
          if (data.flagNorm) for(size_t i=0; i<gridVals.size(); i++) {
              vector<double> bin_sizes;
              bin_sizes.resize(pineappl_grid_bin_count(grid));
              pineappl_grid_bin_normalizations(grid, bin_sizes.data());
              gridVals[i] *= bin_sizes[i];
          }
        }
        // insert values from this grid into output array
        copy_n(gridVals.begin(), gridVals.size(), &val[pos]);
        pos += pineappl_grid_bin_count(grid);
    }
} //compute
