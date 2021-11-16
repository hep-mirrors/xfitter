/*
  @file ReactionAPFELgrid.cc
  @date 2021-11-16
*/

#include"ReactionAPFELgrid.h"
#include"APFELgridGeneration.h"
#include"TFile.h"
#include"xfitter_pars.h"
#include"xfitter_steer.h"
#include"xfitter_cpp_base.h"
#include"appl_grid/appl_grid.h"

#include "APFELgrid/fastkernel.h"                                                                                                                                                                            
#include "APFELgrid/transform.h"                                                                                                                                                                             

#include<memory>
#include"BaseEvolution.h"

using namespace std;
using namespace xfitter;

struct GridData{
  unique_ptr<appl::grid>grid=nullptr;
  int Ndummysize=0;
  //Each grid is either a real APFELgrid or a dummy
  //For real  grids, grid is a valid pointer, and Ndummysize==0
  //For dummy grids, grid==nullptr          , and Ndummysize>0 size of grid
  //Dummy grids return a vector of zeroes of size Ndummysize
};

struct DatasetData{
  vector<GridData>grids;
  int order;
  const double*muR,*muF; // !> renormalisation and factorisation scales
  bool flagNorm;  // !> if true, multiply by bin width
};

const double ONE=1;

extern "C" ReactionAPFELgrid* create() {
  return new ReactionAPFELgrid();
}

void ReactionAPFELgrid::initTerm(TermData*td){
  DatasetData*data=new DatasetData;
  td->reactionData=(void*)data;
  //Split entries in GridName and load grids
  string GridName=td->getParamS("GridName");
  try{
    std::istringstream ss(GridName);
    std::string token;
    while(std::getline(ss, token, ',')){
      data->grids.push_back(GridData());
      GridData&gd=data->grids.back();
      //      appl::grid*g=new appl::grid(token);
    }
  }
  catch ( const std::exception& e ) {
    cerr<<"[FATAL] Unhandled exception while trying to read APFELgrid file(s) \""<<GridName<<"\"; rethrowing exception"<<endl;
    throw e;
  }
  // Get Order
  data->order=OrderMap(td->getParamS("Order"));
  // Get MuR and MuF, or use 1.0 as default
  if(td->hasParam("muR"))data->muR=td->getParamD("muR");
  else data->muR=&ONE;
  if(td->hasParam("muF"))data->muF=td->getParamD("muF");
  else data->muF=&ONE;
  // Get if should normalize by dividing by bin width (no by default)
  data->flagNorm=false;
  if(td->hasParam("norm")){
    int norm=td->getParamI("norm");
    if(norm==1){
      data->flagNorm=true;
    }else if(norm!=0){
      hf_errlog(17102102, "F: unrecognised norm = " + norm);
    }
  }
  size_t Ngrids=data->grids.size();
}



void ReactionAPFELgrid::freeTerm(TermData*td){
  DatasetData*data=(DatasetData*)td->reactionData;
  size_t Ngrids=data->grids.size();
  delete data;
}

void ReactionAPFELgrid::compute(TermData*td,valarray<double>&val,map<string,valarray<double> >&err){
  const DatasetData&data=*(DatasetData*)td->reactionData;
  const int order=data.order;
  const double muR=*data.muR;
  const double muF=*data.muF;
  unsigned int pos = 0;
  //calculate output array size
  {
  size_t np=0;
  for(const GridData&gd:data.grids){
    if(gd.grid)np+=gd.grid->Nobs();
    else       np+=gd.Ndummysize;
  }
  val.resize(np);
  }
  for(const GridData&gd:data.grids){
    appl::grid*grid=gd.grid.get();
    vector<double>gridVals;
    if(grid){//real, non-dummy grid
      gridVals.resize(grid->Nobs());
      td->actualizeWrappers();
      // gridVals=grid->vconvolute(pdf_xfxq_wrapper_,pdf_xfxq_wrapper1_,alphas_wrapper_,order-1,muR,muF,eScale);
      if(data.flagNorm)//scale by bin width
        for(size_t i=0; i<gridVals.size(); i++)
          gridVals[i] *= grid->deltaobs(i);
    }
    // insert values from this grid into output array
    copy_n(gridVals.begin(), gridVals.size(), &val[pos]);
    pos += grid->Nobs();
  }
}
