/*
  @file ReactionAPPLgrid.cc
  @date 2017-03-28
  @author  AddReaction.py
  Created by  AddReaction.py on 2017-03-28
*/

#include"ReactionAPPLgrid.h"
#include"TFile.h"
#include"xfitter_pars.h"
#include"xfitter_steer.h"
#include"xfitter_cpp_base.h"
#include"appl_grid/appl_grid.h"
#include<memory>
#include"BaseEvolution.h"
using namespace std;
using namespace xfitter;
struct GridData{
  unique_ptr<appl::grid>grid=nullptr;
  TH1D*reference=nullptr;//used if flagUseReference=true
  double eScale;// !> CMS energy
  int Ndummysize=0;
  //Each grid is either a real APPLgrid or a dummy
  //For real  grids, grid is a valid pointer, and Ndummysize==0
  //For dummy grids, grid==nullptr          , and Ndummysize>0 size of grid
  //Dummy grids return a vector of zeroes of size Ndummysize
};
struct DatasetData{
  vector<GridData>grids;
  int order;
  const double*muR,*muF; // !> renormalisation and factorisation scales
  bool flagNorm;  // !> if true, multiply by bin width
  bool flagUseReference; // !> if true, prediction will be calculated from reference histogram (for tests and grids validation)
};
const double ONE=1;
extern "C" ReactionAPPLgrid* create() {
  return new ReactionAPPLgrid();
}

void ReactionAPPLgrid::initTerm(TermData*td){
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
      if(beginsWith(token,"DUMMY")){//a dummy grid
        // When GridName=DUMMYX where X is number of bins (e.g. GridName=DUMMY12 for 12 empty bins)
        // reaction will return zeroes (TODO: really?)
        // This is sometimes useful
        //WIP
        gd.Ndummysize=atoi(token.c_str() + 5);//TODO: handle errors
      }else{//a real grid
        appl::grid*g=new appl::grid(token);
        g->trim();
        TFile file(token.c_str());
        TH1D*reference=(TH1D*)file.Get("grid/reference");
        if(!reference)//TODO: maybe do not load reference histogram when !flagUseReference?
          hf_errlog(17033000, "W: no reference histogram grid/reference in " + token);
        else
          reference->SetDirectory(0);//detach the reference histogram from the ROOT file
        gd.grid.reset(g);
        gd.reference=reference;
      }
    }
  }
  catch ( const std::exception& e ) {
    cerr<<"[FATAL] Unhandled exception while trying to read APPLgrid file(s) \""<<GridName<<"\"; rethrowing exception"<<endl;
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
      hf_errlog(17102102, "F: unrecognised norm = " + std::to_string(norm));
      //      hf_errlog(17102102, "F: unrecognised norm = " + norm);
    }
  }
  size_t Ngrids=data->grids.size();
  // Get if should use reference histogram to calculate predictions (for tests or grids validation) (no by default)
  data->flagUseReference=false;
  if(td->hasParam("useReference")){
    int flag=td->getParamI("useReference");
    if(flag==1){
      data->flagUseReference=true;
      // Make sure that the reference histogram is available
      for(size_t i=0;i<Ngrids;++i)
        if(!data->grids[i].reference)
          hf_errlog(17033000, "F: no reference histogram is available");
    }else if(flag!=0){
      hf_errlog(17102102, "F: unrecognised useReference = " + to_string(flag));
    }
  }
  // Get CMS energy (by default the one used to create the grid is used)
  if(td->hasParam("energy")){
    if(data->flagUseReference){
      hf_errlog(17110300,"F: can not apply energy scaling when using predictions from reference histogram");
    }
    double energy=*td->getParamD("energy");
    for(size_t i=0;i<Ngrids;++i){
      double eScale=1;
      double eStored=data->grids[i].grid->getCMSScale();
      if(eStored<1e-3)hf_errlog(17110301, "W: can not apply energy scaling because stored getCMSScale = 0");
      else eScale=eStored/energy;
      data->grids[i].eScale=eScale;
    }
  }else{
    for(size_t i=0;i<Ngrids;++i){
      data->grids[i].eScale=1;
    }
  }
}

void ReactionAPPLgrid::freeTerm(TermData*td){
  DatasetData*data=(DatasetData*)td->reactionData;
  size_t Ngrids=data->grids.size();
  for(size_t i=0;i<Ngrids;++i){
    delete data->grids[i].reference;
  }
  delete data;
}

void  ReactionAPPLgrid::atIteration(){
  double Vud = *XFITTER_PARS::getParamD("Vud");
  double Vus = *XFITTER_PARS::getParamD("Vus");
  double Vub = *XFITTER_PARS::getParamD("Vub");
  double Vcd = *XFITTER_PARS::getParamD("Vcd");
  double Vcs = *XFITTER_PARS::getParamD("Vcs");
  double Vcb = *XFITTER_PARS::getParamD("Vcb");
  double Vtd = *XFITTER_PARS::getParamD("Vtd");
  double Vts = *XFITTER_PARS::getParamD("Vts");
  double Vtb = *XFITTER_PARS::getParamD("Vtb");
  _ckm[0] = { Vud, Vus, Vub };
  _ckm[1] = { Vcd, Vcs, Vcb };
  _ckm[2] = { Vtd, Vts, Vtb };
}

void ReactionAPPLgrid::compute(TermData*td,valarray<double>&val,map<string,valarray<double> >&err){
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
      double eScale=gd.eScale;
      gridVals.resize(grid->Nobs());
      if(!data.flagUseReference){
	auto const grid_ckm =  grid->getckm();
	if (grid_ckm.size()>0) // the grid is for W
	  {
	    // When setting CKM, keep ellements that are zero in the grid zeros (important for V_{tx})
	    std::vector< vector<double> > loc_ckm{ {1,0,0}, {0.,1.,0}, {0.,0.,1} };
	    for ( int irow=0; irow<3; irow+=1)
	      for ( int icol=0; icol<3; icol+=1)
		loc_ckm[irow][icol] =  grid_ckm[irow][icol] == 0 ? 0 : _ckm[irow][icol];
	    grid->setckm(loc_ckm);
	  }
        td->actualizeWrappers();
        gridVals=grid->vconvolute(pdf_xfxq_wrapper_,pdf_xfxq_wrapper1_,alphas_wrapper_,order-1,muR,muF,eScale);
      }else{
        // use reference histogram
        TH1D*ref=gd.reference;
        for(size_t i=0;i<gridVals.size();i++)
          gridVals[i]=ref->GetBinContent(i + 1);
      }
      if(data.flagNorm)//scale by bin width
        for(size_t i=0; i<gridVals.size(); i++)
          gridVals[i] *= grid->deltaobs(i);
    }else{//dummy grid
      gridVals=vector<double>(gd.Ndummysize,0);
    }
    // insert values from this grid into output array
    copy_n(gridVals.begin(), gridVals.size(), &val[pos]);
    pos += grid->Nobs();
  }
  // SZ 27.03.2019 val.size()!=pos should be allowed for bin manipulations
  //if(val.size()!=pos){//TODO: number of data points actually doesn't have to match grid size in some cases, so this check should be replaced by something else
  //  cerr<<"[ERROR] ReactionAPPLgrid: Number of data points ("<<val.size()<<") in term (ID="<<td->id<<") does not match total grid size ("<<pos<<")";
  //  hf_errlog(18072311,"F: ReactionAPPLgrid: number of data points does not match total grid size, see stderr");
  //}
}
