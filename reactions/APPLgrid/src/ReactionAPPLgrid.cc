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
using namespace std;
using namespace xfitter;
 // the class factories
extern "C" ReactionAPPLgrid* create() {
  return new ReactionAPPLgrid();
}
function<void(double const&x,double const&Q,double*pdfs)>active_xfxQ_functions[2];
void xfxWrapper0(const double&x,const double&Q,double*results){active_xfxQ_functions[0](x,Q,results);}
void xfxWrapper1(const double&x,const double&Q,double*results){active_xfxQ_functions[1](x,Q,results);}
ReactionAPPLgrid::ReactionAPPLgrid(){}
ReactionAPPLgrid::~ReactionAPPLgrid(){}
 // Initialize at the start of the computation
int ReactionAPPLgrid::initAtStart(const string &s){return 0;}

 // Initialize for a given dataset:
void ReactionAPPLgrid::setDatasetParameters(int dataSetID, map<string,string> pars, map<string, double> parsDataset) {
  /*DEBUG
  cerr<<"[DEBUG]ReactionAPPLgrid::setDatasetParameters(dataSetID="<<dataSetID<<",pars,parsDataset):\n"
      <<"pars={\n";
  for(map<string,string>::const_iterator it=pars.begin();it!=pars.end();++it){
    cerr<<"  \""<<it->first<<"\":\""<<it->second<<"\"\n";
  }
  cerr<<"}"<<endl;
  cerr<<"parsDataset={\n";
  for(map<string,double>::const_iterator it=parsDataset.begin();it!=parsDataset.end();++it){
    cerr<<"  \""<<it->first<<"\":"<<it->second<<'\n';
  }
  cerr<<"}"<<endl;
  DEBUG*/
  DatasetData&data=dataset_data[dataSetID];
  vector<TH1D*>&references=data.references;
  BaseEvolution*(&evolutions)[2]=data.evolutions;
  map<string,string>::iterator it;
  // Get grid name:
  it=pars.find("GridName");
  if(it!=pars.end()){
    try{
      std::istringstream ss(it->second);
      std::string token;
      while(std::getline(ss, token, ',')){
        std::shared_ptr<appl::grid> g(new appl::grid(token));
        g->trim();
        data.grids.push_back(g);
        TFile file(token.c_str());
        references.push_back((TH1D*)file.Get("grid/reference"));
        if(!references.back())
          hf_errlog(17033000, "W: no reference histogram grid/reference in " + token);
        else
          references.back()->SetDirectory(0);
      }
    }
    catch ( const std::exception& e ) {
      std::string text = "F:Failed to read APPLgrid file(s) "+pars["GridName"];
      hf_errlog_(17032802,text.c_str(),text.size());
    }
  }
  else {
    std::string text = "F:GridName must be specified for the Reaction APPLgrid";
    hf_errlog_(17032801,text.c_str(),text.size());
  }
// Get pointers to evolutions
  for(int i=0;i<2;++i){
    string evolutionName="";
    it=pars.find("evolution"+to_string(i+1));
    if(it!=pars.end())evolutionName=it->second;//else default evolution name will be used
    //cerr<<"[DEBUG]APPLgrid: use evolution \""<<evolutionName<<'\"'<<endl;
    evolutions[i]=get_evolution(evolutionName);
  }
// Determine order
  int order = OrderMap( GetParamS("Order"));  // Global order
  if (pars.find("Order") != pars.end() ) { // Local order
      int localOrder = OrderMap( pars["Order"] );
      order = localOrder>order ? order : localOrder;
  }
  data.order=order;
// Determine MuR and MuF.  Use default 
  data.muR=pars.find("muR") == pars.end() ? GetParam("muR") : stod(pars["muR"]);
  data.muF=pars.find("muF") == pars.end() ? GetParam("muF") : stod(pars["muF"]);

  if(data.muR==0)data.muR=1.0; 
  if(data.muF==0)data.muF=1.0;
  // bin width normalisation (by default no rescaling)
  data.flagNorm=false;
  it = pars.find("norm");
  if (it != pars.end() )
  {
    if(stoi(it->second) == 0)
      data.flagNorm=false;
    else if(stoi(it->second) == 1)
      data.flagNorm=true;
    else
      hf_errlog(17102102, "F: unrecognised norm = " + it->second);
  }
  // use reference histogram to calculate predictions (for tests or grids validation)
  it = pars.find("useReference");
  if (it != pars.end() )
  {
    if(stoi(it->second) == 0)
      data.flagUseReference=false;
    else if(stoi(it->second) == 1)
    {
      data.flagUseReference=true;
      // check that reference histogram is available
      size_t endi=data.references.size();
      for(size_t i=0;i<endi;++i)
        if(!data.references[i])
          hf_errlog(17033000, "W: no reference histogram is available");
    }
    else
      hf_errlog(17102102, "F: unrecognised useReference = " + it->second);
  }
  // CMS energy (by default the one used to create the grid is used)
  it = pars.find("energy");
  unsigned int endi=data.grids.size();
  for(unsigned int i=0;i<endi;++i)
  {
    double eScale = 1.0;
    if (it != pars.end())
    {
      if(data.flagUseReference)
        hf_errlog(17110300, "W: can not apply energy scaling when using predictions from reference histogram");
      else
      {
        double eStored=data.grids[i]->getCMSScale();
        if(eStored < 1e-3)
          hf_errlog(17110301, "W: can not apply energy scaling because stored getCMSScale = 0");
        else
          eScale = eStored/ stof(it->second);
      }
    }
    data.eScale.push_back(eScale);
  }
}


 // Main function to compute results at an iteration
int ReactionAPPLgrid::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err) {
  // iterate over grids
  DatasetData&data=dataset_data.at(dataSetID);
  const int order=data.order;
  const double muR=data.muR;
  const double muF=data.muF;
  BaseEvolution*(&evolutions)[2]=data.evolutions;
  unsigned int pos = 0;
  for(unsigned int g=0;g<data.grids.size();g++)
  {
    auto grid=data.grids[g];
    double eScale=data.eScale[g];
    std::vector<double> gridVals(grid->Nobs());
    if(!data.flagUseReference){
      //For some reason we do not take alphaS from evolutions?
      active_xfxQ_functions[0]=evolutions[0]->xfxQArray();
      if(evolutions[0]==evolutions[1]){
        gridVals=grid->vconvolute(xfxWrapper0,getAlphaS(),order-1,muR,muF,eScale);
      }else{
        active_xfxQ_functions[1]=evolutions[1]->xfxQArray();
        gridVals=grid->vconvolute(xfxWrapper0,xfxWrapper1,getAlphaS(),order-1,muR,muF,eScale);
      }
    }
    else
    {
      // use reference histogram
      for(std::size_t i=0; i<gridVals.size(); i++)
        gridVals[i]=data.references[g]->GetBinContent(i + 1);
    }
    // scale by bin width if requested
    if(data.flagNorm)
      for (std::size_t i=0; i<gridVals.size(); i++)
        gridVals[i] *= grid->deltaobs(i);

    
    // insert values from this grid into output array
    //val.resize(val.size() + grid->Nobs());
    std::copy_n(gridVals.begin(), gridVals.size(), &val[pos]);
    pos += grid->Nobs();
  }
  if(val.size()!=pos){//TODO: number of data points actually doesn't have to match grid size in some cases, so this check should be replaced by something else
    std::ostringstream s;
    s<<"F: ReactionAPPLgrid: Number of data points ("<<val.size()<<") in dataset (ID="<<dataSetID<<") does not match total grid size ("<<pos<<")";
    hf_errlog(18072311,s.str().c_str());
  }
  return 0;
}
