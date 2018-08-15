/*
  @file ReactionAPPLgrid.cc
  @date 2017-03-28
  @author  AddReaction.py
  Created by  AddReaction.py on 2017-03-28
*/

#include "ReactionAPPLgrid.h"
#include "TFile.h"

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
void ReactionAPPLgrid::setDatasetParameters(int dataSetID, map<string,string> pars, map<string, double> parsDataset) {
// Get grid name:
   if ( pars.find("GridName") != pars.end() )  {
     try {
       std::istringstream ss(pars["GridName"]);
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
  int pos = 0;
  for(unsigned int g = 0; g < _grids[dataSetID].size(); g++)
  {
    auto grid = _grids[dataSetID][g];
    std::vector<double> gridVals(grid->Nobs());
    if(!_flagUseReference[dataSetID])
    {
      // Convolute the grid:
      switch (_collType[dataSetID])
      {
        case collision::pp :
          gridVals =  grid->vconvolute( getXFX(), getAlphaS(), _order[dataSetID]-1, _muR[dataSetID], _muF[dataSetID], _eScale[dataSetID][g] );
          break;
        case collision::ppbar :
          gridVals =  grid->vconvolute( getXFX(), getXFX("pbar"), getAlphaS(), _order[dataSetID]-1, _muR[dataSetID], _muF[dataSetID], _eScale[dataSetID][g] );
          break;
        case collision::pn :
          gridVals =  grid->vconvolute( getXFX(), getXFX("n"), getAlphaS(), _order[dataSetID]-1, _muR[dataSetID], _muF[dataSetID], _eScale[dataSetID][g] );
          break;
      }
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
  return 0;
}
