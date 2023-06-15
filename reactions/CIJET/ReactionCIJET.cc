// DB 08/2017 
/*
   @file ReactionCIJET.cc
   @date 2016-12-06
   @author  AddReaction.py
   Created by  AddReaction.py on 2016-12-06
*/

#include "ReactionCIJET.h"


//______________________________________________________________________________
// the class factories
extern "C" ReactionCIJET* create() {
    return new ReactionCIJET();
}
//______________________________________________________________________________
// Initialise for a given dataset:
void ReactionCIJET::initTerm(TermData* td) {
 
    int ID = td->id;  // ID=dataSetID -> updated to termdata td
    if ( td->hasParam("Filename") ) {
        // --- Read file and instantiate CIJET
        const std::string filename = td->getParamS("Filename");
        hf_errlog(17086501,"I: Reading CIJET file "+filename);
        ffnlo.insert(std::make_pair(ID,vector<CIJETReaction*>{new CIJETReaction(filename,this)} ));
    }
    else if ( td->hasParam("Filenames") ) {
        const std::string filenames = td->getParamS("Filenames");
        std::stringstream ss(filenames);
        std::string filename;
        while (std::getline(ss, filename, ',')) {
            hf_errlog(17086501,"I: Reading CIJET file "+filename);
            ffnlo[ID].push_back(new CIJETReaction(filename,this));
        }
    }
    else {
        hf_errlog(17086503,"F:No CIJET file specified. Please provide 'Filename' or 'Filenames' to reaction CIJET.");
        return;
    }
 
    for ( CIJETReaction* fnlo : ffnlo[ID] ) {
        // --- Set order of calculation 
        if ( td->hasParam("Order") ) { // Local order 
        hf_errlog(17094510,"W: Ignoring key 'Order' in .dat file. Only global parameter 'Order' is used.");
        }
        string order = td->getParamS("Order");  // Global order
        bool success=true;
        // CIJET default is 'NLO'
        if      (order=="NNLO") fnlo->setorder(1); // switch NLO ON
        else if (order=="NLO" ) fnlo->setorder(1); // switch NLO ON
        else if (order=="LO"  ) fnlo->setorder(0); // switch NLO OFF
        else    hf_errlog(17094502,"E: CIJET. Unrecognized order: "+order);
        if (!success)  hf_errlog(17094503,"W: CIJET. Requested order "+order+" cannot be set.");
     
        // --- Set CI normalization
        double nm=1.0;
        if ( td->hasParam("CInorm") ) nm = *td->getParamD("CInorm");  // Local order
        else hf_errlog(23061301,"W: Normalization factor 'CInorm' unspecified, assuming 1.0.");
        fnlo->setnorm(nm);
     
        // --- Set scale factors
        double cmur=1, cmuf=1;
        if ( td->hasParam("ScaleFacMuR") ) {
         hf_errlog(17094504,"I: Setting CIJET scale factor mu_R: "+td->getParamS("ScaleFacMuR"));
         cmur = *td->getParamD("ScaleFacMuR");
        }
        if ( td->hasParam("ScaleFacMuF") ) {
            hf_errlog(17094505,"I: Setting CIJET scale factor mu_F: "+td->getParamS("ScaleFacMuF"));
            cmuf = *td->getParamD("ScaleFacMuF");
        }
        fnlo->setscales(cmuf,cmur);
     
        // ---Read the table
        fnlo->setinit();
     
    }
}

//______________________________________________________________________________
// Main function to compute results at an iteration
void ReactionCIJET::compute(TermData* td, valarray<double> &val, map<string, valarray<double> > &err)
{
    td->actualizeWrappers();
    // Get relevant parameters:
    int cicase  = td->getParamI("CIcase");
    double lamotev  = *td->getParamD("LamoTeV");
    double cp1  = *td->getParamD("CI1");
    double cp3  = *td->getParamD("CI3");
    double cp5  = *td->getParamD("CI5");
    // Adjust according to type of contact interactions
    // 1-3 constrained fit
    if ( cicase==1 ) {  //pure left or right handed
        cp3=0.0;   //using cp1 as free parameter
        cp5=0.0;
    }
    else if ( cicase==2 ) { //vector case
        cp3=2*cp1;  //using cp1 as free parameter
        cp5=cp1;
    }
    else if ( cicase==3 ) { //axial-vector case
        cp3=-2*cp1;  //using cp1 as free parameter
        cp5=cp1;
    }
    else  {  //default free-N parameters
    }
 
    vector<double> ci;
    ci.push_back(cp1); 
    ci.push_back(cp3); 
    ci.push_back(cp5);
 
    vector<double> valfnlo;
    valfnlo.reserve(val.size());
    int ID = td->id;
    for ( CIJETReaction* fnlo : ffnlo[ID] ) {
        fnlo->setci(lamotev,ci);
        std::vector<double> cs = fnlo->calcxsec();
        for ( std::size_t i=0; i<cs.size(); i++) valfnlo.push_back(cs[i]);
    }
 
    if ( val.size()!=valfnlo.size() )
        hf_errlog(17116501,"F: Size of CIJET cross section array does not match data.");
    for ( std::size_t i=0; i<val.size(); i++) val[i]=(i<valfnlo.size())?valfnlo[i]:0.;
}
