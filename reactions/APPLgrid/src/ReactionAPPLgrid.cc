/*
  @file ReactionAPPLgrid.cc
  @date 2017-03-28
  @author  AddReaction.py
  Created by  AddReaction.py on 2017-03-28
*/

#include "ReactionAPPLgrid.h"
#include "TFile.h"

#ifdef LHAPDF_ENABLED
LHAPDF::PDF*active_pdf;
void xfxLHAPDF_wrapper(const double&x,const double&Q,double*results){
	vector<double>rtn;
	active_pdf->xfxQ(x,Q,rtn);
	std::copy(rtn.begin(),rtn.end(),results);
}
#endif

 // the class factories
extern "C" ReactionAPPLgrid* create() {
  return new ReactionAPPLgrid();
}
ReactionAPPLgrid::ReactionAPPLgrid(){}
ReactionAPPLgrid::~ReactionAPPLgrid(){
#ifdef LHAPDF_ENABLED
	for(auto it=lhapdf_pdf.begin();it!=lhapdf_pdf.end();++it){
		LHAPDF::PDF*p=it->second;
		if(p)delete p;
	}
#endif
}
 // Initialize at the start of the computation
int ReactionAPPLgrid::initAtStart(const string &s )
{
   return 0;
}

 // Initialize for a given dataset:
void ReactionAPPLgrid::setDatasetParameters(int dataSetID, map<string,string> pars, map<string, double> parsDataset) {
	map<string,string>::iterator it;
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

//Initialize for LHAPDF
//To be used when PDF of second particle is static and of LHAPDF
//Parameters are passed via TermInfo in theory expression in datafile:
//*LHAPDF_SetName
//*LHAPDF_MemberID --- defaults to 0
#ifdef LHAPDF_ENABLED
	{
	it=pars.find("LHAPDF_SetName");
	auto itMemberID=pars.find("LHAPDF_MemberID");
	if(it!=pars.end()){
		int member_id=0;
		try{
			if(itMemberID!=pars.end())member_id=std::stoi(itMemberID->second);
		}catch(const std::exception&ex){
			std::ostringstream s;
			s<<"F:Failed to read PDF memberID for LHAPDF set \""<<it->second<<"\" for use with APPLgrid; dataSetID="<<dataSetID;
			s<<"; Exception:"<<ex.what();
			hf_errlog(3081810,s.str().c_str());
		}
		//Make sure to set collision=LHAPDF
		//std::cout<<"DEBUG LHAPDF_SetName="<<it->second<<std::endl;
		LHAPDF::PDF*p=LHAPDF::mkPDF(it->second,member_id);
		if(!p){
			std::ostringstream s;
			s<<"F:LHAPDF failed to load PDF set \""<<it->second<<"\" for use with APPLgrid; dataSetID="<<dataSetID;
			hf_errlog(28071810,s.str().c_str());
		}
		lhapdf_pdf[dataSetID]=p;
	}else{
		if(itMemberID!=pars.end()){
			std::ostringstream s;
			s<<"W: In APPLgrid reaction parameters (dataSetID="<<dataSetID<<") memberID="<<itMemberID->second<<" is given, but set itself is not";
			hf_errlog(3081811,s.str().c_str());
		}
	}
	}
#else
	if(pars.find("LHAPDF_SetName")!=pars.end()||pars.find("LHAPDF_MemberID")!=pars.end()){
		hf_errlog(25071812,"F:LHAPDF parameter specified for APPLgrid reaction, but xFitter has been compiled without LHAPDF. Use --enable-lhapdf at configure to enable.")
	}
#endif

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
  it=pars.find("collision");
  if (it != pars.end() )
  {
    if     (it->second=="pp")    _collType[dataSetID]=collision::pp;
    else if(it->second=="ppbar") _collType[dataSetID]=collision::ppbar;
    else if(it->second=="pn")    _collType[dataSetID]=collision::pn;
		else if(it->second=="LHAPDF"){
			_collType[dataSetID]=collision::LHAPDF;
			if(lhapdf_pdf.find(dataSetID)==lhapdf_pdf.end())hf_errlog(24071810,"W: collision type=LHAPDF but no LHAPDF set was loaded");
		}
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
      _flagUseReferece[dataSetID] = false;
    else if(stoi(it->second) == 1)
    {
      _flagUseReferece[dataSetID] = true;
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
      if(_flagUseReferece[dataSetID])
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
    if(!_flagUseReferece[dataSetID])
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
#ifdef LHAPDF_ENABLED
        case collision::LHAPDF:
					active_pdf=lhapdf_pdf.find(dataSetID)->second;
          gridVals=grid->vconvolute(getXFX(),xfxLHAPDF_wrapper,getAlphaS(),_order[dataSetID]-1, _muR[dataSetID], _muF[dataSetID], _eScale[dataSetID][g] );
          break;
#endif
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
  if(val.size()!=pos){
    std::ostringstream s;
    s<<"F: Number of data points ("<<val.size()<<") in dataset (ID="<<dataSetID<<") does not match total grid size ("<<pos<<")";
    hf_errlog(18072311,s.str().c_str());
  }
  return 0;
}
