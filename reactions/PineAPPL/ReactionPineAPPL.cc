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
#include <sstream>

using namespace std;
using namespace xfitter;

struct DatasetData {
    vector<pineappl_grid*> grids;
    int order;
    int32_t pdg_id;
    const double *muR, *muF; // !> renormalisation and factorisation scales
    bool flagNorm;  // !> if true, multiply by bin width
    int Nbins;
    int Nord;
    int Nlumi;
    vector<bool> ordervec;
    vector<bool> lumivec;
};

const double ONE=1;

extern "C" ReactionPineAPPL* create() {
    return new ReactionPineAPPL();
}

//AUX function to translate input string into a boolean vector
vector<bool> maskParser(string mask) {
    if (mask.find("[")!=string::npos) mask.erase(0,mask.find("[")+1);
    if (mask.find("]")!=string::npos) mask.erase(mask.find("]"));
    stringstream sstream(mask);
    vector<string> vs;
    vector<bool> vb;
    bool itmp;
    while(sstream.good()) {
        string substr;
        getline(sstream, substr, ',');
        vs.push_back(substr);
    }
    for (string s : vs) {
		if (s!="1" || s!="0") {
		    hf_errlog(23010307, "F: if PineAPPL mask given, it must be a comma-separated string of 1s and 0s!");
        }
        sstream.str("");  sstream.clear();
        sstream << s;
        sstream >> itmp;
        vb.push_back(itmp==0 ? false : true);
    }
    return vb;
}

void ReactionPineAPPL::initTerm(TermData*td) {
    DatasetData* data = new DatasetData;
    td->reactionData = (void*)data;
    // Load grids
    string GridName = td->getParamS("GridName");
    try {
        istringstream ss(GridName);
        string token;
        while (getline(ss, token, ',')) {
            pineappl_grid* g = pineappl_grid_read(token.c_str());             
            data->grids.push_back(g);
        }
        // Init dimensions
        data->Nbins = pineappl_grid_bin_count(data->grids[0]);
        data->Nord  = pineappl_grid_order_count(data->grids[0]);
        auto* lumi  = pineappl_grid_lumi(data->grids[0]);
        data->Nlumi = pineappl_lumi_count(lumi);
        hf_errlog(22122901, "I: read PineAPPL grids with "
                            +to_string(data->Nbins)
                            +" bins");
    } catch ( const exception& e ) {
        hf_errlog(22121301, "E: Unhandled exception while trying to read PineAPPL grid(s) "+GridName);
        throw e;
    }

    if (td->hasParam("PDG")) data->pdg_id = td->getParamI("PDG");
    else {
        data->pdg_id = 2212;
        hf_errlog(22121201, "W: PineAPPL initTerm assuming protons");
    }

    // Get OrderMask
    if (td->hasParam("OrderMask")) {
        string ordS = td->getParamS("OrderMask");
        data->ordervec = maskParser(ordS);
        hf_errlog(23010305, "I: PineAPPL order mask: "+ordS);
        if (data->Nord!=data->ordervec.size()) {
            hf_errlog(23010301, "F: PineAPPL grid order count: "
                                +to_string(data->Nord)
                                +" disagrees with order mask size: "
                                +to_string(data->ordervec.size()));
        }
    } else data->Nord = 0;
    if (data->Nord==0) hf_errlog(23010302, "I: PineAPPL order mask unspecified, using all in grid");

	// Get LumiMask
    if (td->hasParam("LumiMask")) {
        string lumiS = td->getParamS("LumiMask");
        data->lumivec = maskParser(lumiS);
        if (data->Nlumi!=data->lumivec.size()) {
            hf_errlog(23010303, "F: PineAPPL grid lumi count: "
                                +to_string(data->Nlumi)
                                +" disagrees with lumi mask size: "
                                +to_string(data->lumivec.size()));
        }
    } else data->Nlumi = 0;
    if (data->Nlumi==0) hf_errlog(23010304, "I: PineAPPL lumi mask unspecified, using all in grid");

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
    const double muR = *data.muR;
    const double muF = *data.muF;
    unsigned int pos = 0;
    bool order_mask[data.Nord];
    bool lumi_mask[data.Nlumi];
    for (int i=0; i<data.Nord; ++i) order_mask[i] = data.ordervec[i];
    for (int i=0; i<data.Nlumi; ++i) lumi_mask[i] = data.lumivec[i];   

    //calculate output array size
    size_t np=0;
    for (pineappl_grid* grid : data.grids) if (grid) np += pineappl_grid_bin_count(grid);
    val.resize(np);

    for (pineappl_grid* grid : data.grids) {
        vector<double> gridVals;
        gridVals.resize(data.Nbins);
        if (grid) {//real, non-dummy grid
            td->actualizeWrappers();
            //See function specification in deps/pineappl/include/pineappl_capi/pineappl_capi.h
            pineappl_grid_convolute_with_one(grid, data.pdg_id, 
                                             &xfx, &alphas, 
                                             nullptr,//"state" provided to wrappers, redundant in xFitter
                                             data.Nord>0 ? order_mask : nullptr,
                                             data.Nlumi>0 ? lumi_mask : nullptr,
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
