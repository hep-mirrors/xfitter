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
ReactionAPPLgrid::ReactionAPPLgrid(){
  evolutions[0]=nullptr;
  evolutions[1]=nullptr;
}
ReactionAPPLgrid::~ReactionAPPLgrid(){}
 // Initialize at the start of the computation
int ReactionAPPLgrid::initAtStart(const string &s){return 0;}

 // Initialize for a given dataset:
void ReactionAPPLgrid::setDatasetParameters(int dataSetID, map<string,string> pars, map<string, double> parsDataset) {
  //<DEBUG
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
  //DEBUG>
	map<string,string>::iterator it;
	// Get grid name:
  it=pars.find("GridName");
  if(it!=pars.end()){
    try{
      std::istringstream ss(it->second);
       //std::cout << pars["GridName"] << '\n';
       std::string token;
       while(std::getline(ss, token, ','))
       {
         //std::cout << token << '\n';
         std::shared_ptr<appl::grid>  g(new appl::grid(token));
         g->trim();
         _grids[dataSetID].push_back(g);
         TFile file(token.c_str());
         _references[dataSetID].push_back((TH1D*)file.Get("grid/reference"));
         if(!_references[dataSetID].back())
           hf_errlog(17033000, "W: no reference histogram grid/reference in " + token);
         else
           _references[dataSetID].back()->SetDirectory(0);
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
    cerr<<"[DEBUG]APPLgrid: use evolution \""<<evolutionName<<'\"'<<endl;
    evolutions[i]=get_evolution(evolutionName);
  }
// Determine order
   int order = OrderMap( GetParamS("Order"));  // Global order
   if (pars.find("Order") != pars.end() ) { // Local order
       int localOrder = OrderMap( pars["Order"] );
       order = localOrder>order ? order : localOrder;
   }

   _order[dataSetID] = order;
// Determine MuR and MuF.  Use default 
  _muR[dataSetID] =  pars.find("muR") == pars.end() ? GetParam("muR") : stod(pars["muR"]);
  _muF[dataSetID] =  pars.find("muF") == pars.end() ? GetParam("muF") : stod(pars["muF"]);

  if (_muR[dataSetID] == 0) _muR[dataSetID] = 1.0; 
  if (_muF[dataSetID] == 0) _muF[dataSetID] = 1.0;
  // bin width normalisation (by default no rescaling)
  _flagNorm[dataSetID] = false;
  it = pars.find("norm");
  if (it != pars.end() )
  {
    if(stoi(it->second) == 0)
      _flagNorm[dataSetID] = false;
    else if(stoi(it->second) == 1)
      _flagNorm[dataSetID] = true;
    else
      hf_errlog(17102102, "F: unrecognised norm = " + it->second);
  }
  // use reference histogram to calculate predictions (for tests or grids validation)
  it = pars.find("useReference");
  if (it != pars.end() )
  {
    if(stoi(it->second) == 0)
      _flagUseReference[dataSetID] = false;
    else if(stoi(it->second) == 1)
    {
      _flagUseReference[dataSetID] = true;
      // check that reference histogram is available
      for(std::size_t i=0; i<_references[dataSetID].size(); i++)
        if(!_references[dataSetID][i])
          hf_errlog(17033000, "W: no reference histogram is available");
    }
    else
      hf_errlog(17102102, "F: unrecognised useReference = " + it->second);
  }
  // CMS energy (by default the one used to create the grid is used)
  it = pars.find("energy");
  for(unsigned int i = 0; i < _grids[dataSetID].size(); i++)
  {
    double eScale = 1.0;
    if (it != pars.end())
    {
      if(_flagUseReference[dataSetID])
        hf_errlog(17110300, "W: can not apply energy scaling when using predictions from reference histogram");
      else
      {
        double eStored = _grids[dataSetID][i]->getCMSScale();
        if(eStored < 1e-3)
          hf_errlog(17110301, "W: can not apply energy scaling because stored getCMSScale = 0");
        else
          eScale = _grids[dataSetID][i]->getCMSScale() / stof(it->second);
      }
    }
    _eScale[dataSetID].push_back(eScale);
  }
}


 // Main function to compute results at an iteration
int ReactionAPPLgrid::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err) {
  // iterate over grids
  unsigned int pos = 0;
  for(unsigned int g = 0; g < _grids[dataSetID].size(); g++)
  {
    auto grid = _grids[dataSetID][g];
    std::vector<double> gridVals(grid->Nobs());
    if(!_flagUseReference[dataSetID]){
      //For some reason we do not take alphaS from evolutions?
      active_xfxQ_functions[0]=evolutions[0]->xfxQArray();
      active_xfxQ_functions[1]=evolutions[1]->xfxQArray();
      gridVals=grid->vconvolute(xfxWrapper0,xfxWrapper1,getAlphaS(),_order[dataSetID]-1,_muR[dataSetID],_muF[dataSetID],_eScale[dataSetID][g]);
      //TODO: call specialized vconvolute when the two evolutions are the same
    }
    else
    {
      // use reference histogram
      for(std::size_t i=0; i<gridVals.size(); i++)
        gridVals[i] = _references[dataSetID][g]->GetBinContent(i + 1);
    }
    // scale by bin width if requested
    if(_flagNorm[dataSetID])
      for (std::size_t i=0; i<gridVals.size(); i++)
        gridVals[i] *= grid->deltaobs(i);

    
    // insert values from this grid into output array
    std::copy_n(gridVals.begin(), gridVals.size(), &val[pos]);
    pos += grid->Nobs();
  }
  if(val.size()!=pos){
    std::ostringstream s;
    s<<"F: Number of data points ("<<val.size()<<") in dataset (ID="<<dataSetID<<") does not match total grid size ("<<pos<<")";
    hf_errlog(18072311,s.str().c_str());
  }
  return 0;
}
