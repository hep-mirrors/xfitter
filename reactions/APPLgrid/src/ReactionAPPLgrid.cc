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

 // Initialize at the start of the computation
int ReactionAPPLgrid::initAtStart(const string &s )
{
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

// Determine if pp or ppbar
  _collType[dataSetID] = collision::pp;
  if (parsDataset.find("ppbar") != parsDataset.end() ) {
    _collType[dataSetID] = collision::ppbar;
  }
  if (parsDataset.find("pn") != parsDataset.end() ) {
    _collType[dataSetID] = collision::pn;
  }
  // check if collision settings are provided in the new format key=value
  map<string,string>::iterator it = pars.find("collision");
  if (it != pars.end() )
  {
    if(it->second == "pp")
      _collType[dataSetID] = collision::pp;
    else if(it->second == "ppbar")
      _collType[dataSetID] = collision::ppbar;
    else if(it->second == "pn")
      _collType[dataSetID] = collision::pn;
    else
      hf_errlog(17102101, "F: unrecognised collision type = " + it->second);
  }
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
}


 // Main function to compute results at an iteration
int ReactionAPPLgrid::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err) {
 // Convolute the grid:
  std::vector<double> vals(val.size());
  switch (_collType[dataSetID])
    {
    case collision::pp :  
      vals =  _grids[dataSetID]->vconvolute( getXFX(), getAlphaS(), _order[dataSetID]-1, _muR[dataSetID], _muF[dataSetID] );
      break; 
    case collision::ppbar :  
      vals =  _grids[dataSetID]->vconvolute( getXFX(), getXFX("pbar"), getAlphaS(), _order[dataSetID]-1, _muR[dataSetID], _muF[dataSetID] );
      break;
    case collision::pn :
      vals =  _grids[dataSetID]->vconvolute( getXFX(), getXFX("n"), getAlphaS(), _order[dataSetID]-1, _muR[dataSetID], _muF[dataSetID] );
      break;
    }
  for (std::size_t i=0; i<vals.size(); i++) {
    val[i] = vals[i];
    // scale by bin width if requested
    if(_flagNorm[dataSetID])
    {
      double bw = _grids[dataSetID]->deltaobs(i + 1);
      val[i] *= bw;
    }
  }
  return 0;
}

