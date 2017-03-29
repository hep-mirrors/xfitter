/*
  @file ReactionAPPLgrid.cc
  @date 2017-03-28
  @author  AddReaction.py
  Created by  AddReaction.py on 2017-03-28
*/

#include "ReactionAPPLgrid.h"

 // the class factories
extern "C" ReactionAPPLgrid* create() {
  return new ReactionAPPLgrid();
}


// Place here temporary. Should go in simplified xfitter_cpp.h-like.
extern "C" {
  void hf_errlog_(const int &id, const char text[], int);
}
int OrderMap(std::string ord) {
  std::map <string,int> mOrderMap  = { {"NNNLO",4}, {"NNLO",3}, {"NLO",2}, {"LO",1} };
  if ( mOrderMap.find(ord) != mOrderMap.end() ) {
    return mOrderMap[ord];
  }
  else {
    string text = "F: Can not convert "+ord+" to computation order";
    hf_errlog_(17032803,text.c_str(),text.size());
    return -1;
  }
}


 // Initialize at the start of the computation
int ReactionAPPLgrid::initAtStart(const string &s )
{
  //   gXFX = getXFX();
  //  gAs  = alpha_S;
   return 0;
}

 // Initialisze for a given dataset:
void ReactionAPPLgrid::setDatasetParamters(int dataSetID, map<string,string> pars, map<string, double> parsDataset) {
// Get grid name:
   if ( pars.find("GridName") != pars.end() )  {
     try {
        std::shared_ptr<appl::grid>  g(new appl::grid(pars["GridName"])); 
        g->trim();
        _grids[dataSetID] = g;
     }
     catch ( const std::exception& e ) {
        std::string text = "F:Failed to read APPLgrid file"+pars["GridName"];
        hf_errlog_(17032802,text.c_str(),text.size());
     }
   }
   else {
     std::string text = "F:GridName must be specified for the Reaction APPLgrid";
     hf_errlog_(17032801,text.c_str(),text.size());
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

}


 // Main function to compute results at an iteration
int ReactionAPPLgrid::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err) {
 // Convolute the grid:
   std::vector<double> vals =  _grids[dataSetID]->vconvolute( getXFX(), getAlphaS(), _order[dataSetID]-1, _muR[dataSetID], _muF[dataSetID] );

   int j=0; 
   for (std::size_t i=0; i<vals.size(); i++) {
      if ( notMasked(dataSetID, i))  {
	 val[j++] = vals[i];
      }
   }
  return 0;
}

