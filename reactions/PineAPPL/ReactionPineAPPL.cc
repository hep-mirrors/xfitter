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
#include <unistd.h>
#include <sys/shm.h>
#include <sys/wait.h>
#include "BaseEvolution.h"
#include "TermData.h"
#include <sstream>
#include <numeric>

using namespace std;
using namespace xfitter;

struct DatasetData {
    vector<pineappl_grid*> grids;
    int order;
    const double *muR, *muF; // !> renormalisation and factorisation scales
    bool flagNorm;  // !> if true, multiply by bin width
    int Nord;
    int Nlumi;
    double energyRescale;
    vector<bool> ordervec;
    vector<bool> lumivec;
    std::string evolution;
    std::string evolution1;
    std::string evolution2;
    vector<std::string> GridNames;
    std::vector<std::vector<int> > rebin;
    const std::string get_pars_as_str() const {
        auto ordervec_rev = ordervec;
        reverse(ordervec_rev.begin(), ordervec_rev.end());
        auto order_key = accumulate(ordervec_rev.rbegin(), ordervec_rev.rend(), 0, [](int x, int y) { return (x << 1) + y; });
        auto lumivec_rev = ordervec;
        reverse(lumivec_rev.begin(), lumivec_rev.end());
        auto lumi_key = accumulate(lumivec_rev.rbegin(), lumivec_rev.rend(), 0, [](int x, int y) { return (x << 1) + y; });
        const std::string key = std::to_string(energyRescale) + "_" + std::to_string((unsigned long long)(void**)muR) + "_" + std::to_string((unsigned long long)(void**)muF) +
            "_" + std::to_string(Nord) + "_" + std::to_string(order_key) + "_" + std::to_string(Nlumi) + "_" + std::to_string(lumi_key) + 
            "_" + std::to_string(flagNorm) + "_" + evolution + "_" + evolution1 + "_" + evolution2;
        return key;
    }
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
        if (s!="1" && s!="0") {
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
    super::initTerm(td);
    DatasetData* data = new DatasetData;
    td->reactionData = (void*)data;
    // Load grids
    string GridName = td->getParamS("GridName");
    try {
        istringstream ss(GridName);
        string token;
        while (getline(ss, token, ',')) {
            pineappl_grid* g = nullptr;
            const auto& find = _initialized.find(token);
            if(find == _initialized.end()) {
                g = pineappl_grid_read(token.c_str());             
                _initialized.insert(std::make_pair(token, g));
                hf_errlog(22122901, "I: read PineAPPL grids with "
                                    +to_string(pineappl_grid_bin_count(g))
                                    +" bins");
            }
            else {
                g = find->second;
            }
            data->grids.push_back(g);
            data->GridNames.push_back(token);
        }
        // Init dimensions
        data->Nord  = pineappl_grid_order_count(data->grids[0]);
        auto* lumi  = pineappl_grid_lumi(data->grids[0]);
        data->Nlumi = pineappl_lumi_count(lumi);
    } catch ( const exception& e ) {
        hf_errlog(22121301, "E: Unhandled exception while trying to read PineAPPL grid(s) "+GridName);
        throw e;
    }

    // Get OrderMask
    if (td->hasParam("OrderMask")) {
        string ordS = td->getParamS("OrderMask");
        data->ordervec = maskParser(ordS);
        hf_errlog(23010305, "I: PineAPPL order mask: "+ordS);
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

    // Get CMS energy (by default the one used to create the grid is used)
    if(td->hasParam("energyRescale")) {
        data->energyRescale = *td->getParamD("energyRescale");
    }
    else {
        data->energyRescale = 1.0;
    }

    // evolution
    if(td->hasParam("evolution")) 
        data->evolution = td->getParamS("evolution");
    else
        data->evolution = "";
    if(td->hasParam("evolution1"))
        data->evolution1 = td->getParamS("evolution1");
    else
        data->evolution1 = "";
    if(td->hasParam("evolution2"))
        data->evolution2 = td->getParamS("evolution2");
    else
        data->evolution2 = "";

    // rebin
    if (td->hasParam("rebin")) {
        std::vector<std::string> rebin_vars;
        std::vector<double> rebin_vars_bound;
        string rebin_str = td->getParamS("rebin");
        istringstream ss(rebin_str);
        string token;
        while (getline(ss, token, ',')) {
            auto start = token.find("[");
            if (start != std::string::npos) {
                rebin_vars.push_back(std::string(token.data(), start));
                auto start2 = token.find(",");
                if (token[token.length()-1] != ']')
                    hf_errlog(23060601, "F: wrong format of rebin: " + rebin_str);
                auto bound_str = std::string(token.data() + start + 1);
                bound_str = std::string(bound_str.data(), bound_str.length()-1);
                rebin_vars_bound.push_back(atof(bound_str.data()));
            }
            else
                rebin_vars.push_back(token);
                rebin_vars_bound.push_back(std::numeric_limits<double>::quiet_NaN());
        }
        if (data->grids.size() == 0)
            hf_errlog(23060501, "F: cannot rebin without grids");
        size_t np = 0;
        int ndim = -1;
        for (size_t ig = 0; ig < data->grids.size(); ig++) {
            if (data->grids[ig]) {
                np += pineappl_grid_bin_count(data->grids[ig]);
                auto thisdim = pineappl_grid_bin_dimensions(data->grids[ig]);
                if (ndim == -1)
                ndim = thisdim;
                else if (ndim != thisdim)
                hf_errlog(23060502, "F: Dimension mismatch for grid " + data->GridNames[ig]);
            }
        }
        if (rebin_vars.size() != (ndim * 2))
            hf_errlog(23060502, "F: rebin [" + std::to_string(rebin_vars.size()) + "] and grid dimension [" + std::to_string(ndim) + "] mismatch");
        std::vector<std::vector<double> > binsl(ndim);
        std::vector<std::vector<double> > binsr(ndim);
        for (int idim = 0; idim < ndim; idim++){
            binsl[idim].resize(np);
            binsr[idim].resize(np);
        }
        unsigned int pos = 0;
        for (size_t igrid = 0; igrid < data->grids.size(); igrid++) {
            pineappl_grid* grid = data->grids[igrid];
            if (grid) {//real, non-dummy grid
                for (int idim = 0; idim < ndim; idim++){
                    pineappl_grid_bin_limits_left(grid, idim, binsl[idim].data() + pos);
                    pineappl_grid_bin_limits_right(grid, idim, binsr[idim].data() + pos);
                }
                pos += pineappl_grid_bin_count(grid);
            }
        }
        std::vector<std::vector<double> > bins_data(rebin_vars.size());
        for (size_t ivar = 0; ivar < rebin_vars.size(); ivar++) {
            const auto& thisbins = *(const_cast<std::valarray<double>*>(td->getBinColumnOrNull(rebin_vars[ivar])));
            bins_data[ivar].resize(thisbins.size());
            std::copy(&thisbins[0], &thisbins[0] + thisbins.size(), bins_data[ivar].begin());
        }
        auto compare = [](double b1, double b2) {
            const double eps = 1e-6;
            auto diff = b2 - b1;
            auto reldiff = (b1 == 0.) ? ( (b2 == 0.) ? 0. : 1. ) : diff / b1;
            if (fabs(diff) < eps || fabs(reldiff) < eps) return 0;
            else if (diff > 0) return 1;
            else return -1;
        };
        auto& rebin = data->rebin;
        rebin.resize(bins_data[0].size());
        for(size_t bin = 0; bin < bins_data[0].size(); bin++) {
            for (int idim = 0; idim < ndim; idim++) {
                if (!std::isnan(rebin_vars_bound[0+idim*2]) && bins_data[0+idim*2][bin] <= rebin_vars_bound[0+idim*2])
                    bins_data[0+idim*2][bin] = binsl[idim][0];
                if (!std::isnan(rebin_vars_bound[1+idim*2]) && bins_data[1+idim*2][bin] <= rebin_vars_bound[1+idim*2])
                    bins_data[1+idim*2][bin] = binsr[idim][binsr[idim].size()-1];
            }
            rebin[bin].resize(binsl[0].size());
            std::vector<int> match_l(ndim, 0);
            std::vector<int> match_r(ndim, 0);
            for(size_t bingrid = 0; bingrid < binsl[0].size(); bingrid++) {
                rebin[bin][bingrid] = 0;
                std::vector<int> flag(ndim);
                for (size_t idim = 0; idim < ndim; idim++) {
                    flag[idim] = 0;
                    auto l = compare(bins_data[0+idim*2][bin], binsl[idim][bingrid]);
                    if (l == -1) continue;
                    else {
                        if (l == 0) match_l[idim] = 1;
                        auto r = compare(bins_data[1+idim*2][bin], binsr[idim][bingrid]);
                        if (r == 1) continue;
                        if (r == 0) match_r[idim] = 1;
                        flag[idim] = 1;
                    }
                }
                if (std::all_of(flag.begin(), flag.end(), [](bool v) { return v; }))
                    rebin[bin][bingrid] = 1;
            }
            for (size_t idim = 0; idim < ndim; idim++) {
                if (match_l[idim] == 0) hf_errlog(23051203, "F: Binning mismatch for " + rebin_vars[0+idim*2] + " " + std::to_string(bins_data[0+idim*2][bin]));
                if (match_r[idim] == 0) hf_errlog(23051203, "F: Binning mismatch for " + rebin_vars[1+idim*2] + " " + std::to_string(bins_data[1+idim*2][bin]));
            }
        }
    }

    // store all grids to convolute them in atIteration
    const std::string key_pars = data->get_pars_as_str();
    for (const auto&  g : data->grids) {
        std::string key = std::to_string((unsigned long long)(void**)g) + "_" + key_pars;
        if (_convolved.find(key) == _convolved.end()) {
            _convolved[key] = std::make_pair(vector<double>(), td);
            _convolved_vector_of_keys.push_back(key);
        }
    }
} //initTerm

void ReactionPineAPPL::freeTerm(TermData*td) {
    DatasetData* data = (DatasetData*)td->reactionData;
    size_t Ngrids=data->grids.size();
    for (size_t i=0; i<Ngrids; ++i) pineappl_grid_delete(data->grids[i]);
    delete data;
}

void ReactionPineAPPL::atIteration() {
    //Pineappl assumes PDF and alphaS wrapper function pointers
    //are given a pointer "state", e.g. an LHAPDF object, if the 
    //values were read from there. Here, call existing wrappers 
    //rearranging parameter list & return value to suit pineappl 
    //convolution function.
    //printf("atIteration _pacount = %d\n", _pacount);
    _pacount = 0;
    auto xfx = [](int32_t id_in, double x, double q2, void *state) {
        if(x >= 1.) {
            return 0.;
        }
        double pdfs[13];
        int32_t id = id_in==21 ? 6 : id_in+6;
        double energyRescale = *((double*)state);
        auto x_actual = x*energyRescale;
        if(x_actual >= 1.) {
            return 0.;
        }
        pdf_xfxq_wrapper_(x_actual, sqrt(q2), pdfs);
        return pdfs[id];
    };
    auto xfx1 = [](int32_t id_in, double x, double q2, void *state) {
        double pdfs[13];
        int32_t id = id_in==21 ? 6 : id_in+6;
        double energyRescale = *((double*)state);
        auto x_actual = x*energyRescale;
        if(x_actual >= 1.) {
            return 0.;
        }
        pdf_xfxq_wrapper1_(x_actual, sqrt(q2), pdfs);
        return pdfs[id];
    };
    auto alphas = [](double q2, void *state) {
        return alphas_wrapper_(sqrt(q2));
    };
    // Fix PDG ID to p, avoiding double charge conjugation in case pbar is used,
    // as this is already done elsewhere before passing PDFs to PineAPPL
    int32_t PDGID = 2212;  //DO NOT MODIFY

    auto calc_one = [&](int i, std::vector<double>& gridVals) {
        const std::string& key = _convolved_vector_of_keys[i];
        auto& it = _convolved[_convolved_vector_of_keys[i]];
        TermData* td = it.second;
        td->actualizeWrappers();
        const DatasetData& data = *(DatasetData*)td->reactionData;
        bool order_mask[data.Nord];
        bool lumi_mask[data.Nlumi];
        for (int i=0; i<data.Nord; ++i) order_mask[i] = data.ordervec[i];
        for (int i=0; i<data.Nlumi; ++i) lumi_mask[i] = data.lumivec[i];
        pineappl_grid* grid = (pineappl_grid*)stoull(key.substr(0, key.std::string::find('_')));
        gridVals.resize(pineappl_grid_bin_count(grid));
        //See function specification in deps/pineappl/include/pineappl_capi/pineappl_capi.h
        _pacount++;
        pineappl_grid_convolute_with_two(grid, 
                                            PDGID, xfx, 
                                            PDGID, xfx1, 
                                            alphas, 
                                            (void*)&data.energyRescale,//"state" is energyRescale
                                            data.Nord>0 ? order_mask : nullptr,
                                            data.Nlumi>0 ? lumi_mask : nullptr,
                                            *data.muR, *data.muF, gridVals.data());
        //scale by bin width
        if (data.flagNorm) for(size_t i=0; i<gridVals.size(); i++) {
            vector<double> bin_sizes;
            bin_sizes.resize(pineappl_grid_bin_count(grid));
            pineappl_grid_bin_normalizations(grid, bin_sizes.data());
            gridVals[i] *= bin_sizes[i]*std::pow(data.energyRescale, 2.);
        }
    };

    int ngrids = _convolved.size();
    int ncpu =  xfitter::xf_ncpu(_ncpu);
    if (ncpu == 1) {
        for (int i = 0; i < ngrids; i++) {
            calc_one(i, _convolved[_convolved_vector_of_keys[i]].first);
        }
    }
    else {
        std::vector<int> positions(ngrids);
        int nbins = 0;
        for (size_t igrid = 0; igrid < _convolved_vector_of_keys.size(); igrid++) {
            positions[igrid] = nbins;
            const std::string& key = _convolved_vector_of_keys[igrid];
            pineappl_grid* grid = (pineappl_grid*)stoull(key.substr(0, key.std::string::find('_')));
            nbins += pineappl_grid_bin_count(grid);
        }
        // Shared memory for predictions
        int shmid;
        double* sharedArray;
        shmid = shmget(IPC_PRIVATE, sizeof(double) * nbins, IPC_CREAT | 0666);
        if (shmid < 0) {
        hf_errlog(2023060200,"F: Failed to create shared memory segment");
        }
        sharedArray = static_cast<double*>(shmat(shmid, nullptr, 0));
        if (sharedArray == reinterpret_cast<double*>(-1)) {
        hf_errlog(2023060201,"F: Failed to attach shared memory segment");
        }
        // define Chunks
        int chunkSize = ngrids / ncpu;
        int reminder  = ngrids % ncpu; 
        int first = 0;
        int startIndex = 0;
        int endIndex = 0;
        // loop over all
        for (int icpu = 0; icpu < min(ncpu, ngrids); icpu++) {
            startIndex = endIndex;
            endIndex   = startIndex + chunkSize;
            if (icpu < reminder) {
                endIndex += 1;
            }
            pid_t pid = xfitter::xf_fork( min(ncpu, ngrids)  );
            if ( pid == 0) {       
                    // close all open files (e.g. minuit.out.txt) to avoid multiple buffered output
                    int fdlimit = (int)sysconf(_SC_OPEN_MAX);
                    for (int i = STDERR_FILENO + 1; i < fdlimit; i++) {
                    close(i);
                }
                for (int i = first+startIndex; i < first+endIndex; i++) {
                    //printf("CPU %d computing %d\n", icpu, i);
                    std::vector<double> gridVals = std::vector<double>();
                    calc_one(i, gridVals);
                    //printf("computed\n");
                    const std::string& key = _convolved_vector_of_keys[i];
                    pineappl_grid* grid = (pineappl_grid*)stoull(key.substr(0, key.std::string::find('_')));
                    for (size_t ibin = 0; ibin < gridVals.size(); ibin++)
                        sharedArray[positions[i] + ibin] = gridVals[ibin];
                }
                exit(0);	    
            }
            else if (pid<0) {
                hf_errlog(2023060204,"F: Failed to create a fork process");	
            }
        }	
        // Wait ...
        int status;
        while (wait(&status) > 0);    
        // Store result
        for (size_t i = 0; i<ngrids; i++) {
            const std::string& key = _convolved_vector_of_keys[i];
            pineappl_grid* grid = (pineappl_grid*)stoull(key.substr(0, key.std::string::find('_')));
            nbins = pineappl_grid_bin_count(grid);
            _convolved[_convolved_vector_of_keys[i]].first.resize(nbins);
            for(int ibin = 0; ibin < nbins; ibin++) {
                _convolved[_convolved_vector_of_keys[i]].first[ibin] = sharedArray[positions[i] + ibin];
            }
        }    
        // Detach and remove shared memory segments
        shmdt(sharedArray);
        shmctl(shmid, IPC_RMID, NULL);
    }
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
    //err.resize(np);

    const std::string key_pars = data.get_pars_as_str();
    for (size_t igrid = 0; igrid < data.grids.size(); igrid++) {
        pineappl_grid* grid = data.grids[igrid];
        std::string key = std::to_string((unsigned long long)(void**)grid) + "_" + key_pars;
        vector<double> gridVals = _convolved[key].first;
        // insert values from this grid into output array
        copy_n(gridVals.begin(), gridVals.size(), &val[pos]);
        pos += pineappl_grid_bin_count(grid);
    }
    // rebin
    if (data.rebin.size() > 0) {
        const auto& rebin = data.rebin;
        auto val_orig = val;
        val.resize(rebin.size());
        //auto err_orig = err;
        //err.resize(rebin.size());
        for(size_t i1 = 0; i1 < rebin.size(); i1++) {
            val[i1] = 0.;
            //err[i1] = 0.;
            for (size_t i2 = 0; i2 < rebin[i1].size(); i2++) {
                val[i1] += val_orig[i2] * rebin[i1][i2];
                //err[i1] += err_orig[i2] * rebin[i1][i2];
            }
        }
    }
} //compute
