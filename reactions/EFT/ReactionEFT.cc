   /*
     @file ReactionEFT.cc
     @date 2023-05
   */

#include "ReactionEFT.h"

//______________________________________________________________________________
// the class factories
extern "C" ReactionEFT* create() {
  return new ReactionEFT();
}

//______________________________________________________________________________
// Initialise for a given dataset:
void ReactionEFT::initTerm(TermData* td) {
  vector<string> fname_list;

  int ID = td->id;  // ID=dataSetID -> updated to termdata td
  int find_input_Q = 0;

  //------------------------------------------------------------------
  // construct EFTReaders
  if ( td->hasParam("Filename") ) {
    find_input_Q = 1;

    const std::string filename = td->getParamS("Filename");
    hf_errlog(23032801,"I: Reading EFT file "+filename);
    fname_list.push_back(filename);
  }

  if ( td->hasParam("Filenames") ) {
    if (find_input_Q > 0) 
      hf_errlog(23052301, "F: there can be only one input in Filename/Filenames/MixedInput");
    find_input_Q = 1;

    const std::string filenames = td->getParamS("Filenames");
    std::stringstream ss(filenames);
    std::string filename;
    while (std::getline(ss, filename, ',')) {
      hf_errlog(23032801,"I: Reading EFT file "+filename);
      fname_list.push_back(filename);
    }
  }

  if ( td->hasParam("MixedInput") ) {
    if (find_input_Q > 0) 
      hf_errlog(23052301, "F: there can be only one input in Filename/Filenames/MixedInput");
    find_input_Q = 2;

    const std::string filename = td->getParamS("MixedInput");
    hf_errlog(23032801,"I: Reading EFT file "+filename);
    fname_list.push_back(filename);
  }

  if (find_input_Q == 0)
    hf_errlog(23032803,"F:No EFT file specified. Input one of 'Filename/Filenames/MixedInput' for EFT terms.");
  else if (find_input_Q == 1)
    EFT_terms.insert(std::make_pair(ID, new EFTReader(fname_list, "fixed", this) ));
  else
    EFT_terms.insert(std::make_pair(ID, new EFTReader(fname_list, "mixed", this) ));


  //------------------------------------------------------------------
  // read the list of EFT parameters
  // todo: what if ListEFTParam is not available
  string list_EFT_param = td->getParamS("ListEFTParam");
  vector<string> name_EFT_param;

  if (list_EFT_param.size() == 0) {
    hf_errlog(23040601,"I: list of EFT parameters is empty");
  } else {
    stringstream ss(list_EFT_param);
    string param_name;
    while(getline(ss, param_name, ',')) {
      name_EFT_param.push_back(param_name);
    }
  }

  EFT_terms[ID]->initParamName(name_EFT_param);

  //-------------------------------------------------------
  // more init
  if ( td->hasParam("Debug") )
    EFT_terms[ID]->debug = td->getParamI("Debug");

  if ( td->hasParam("AbsOutput") ) {
    string s1 = td->getParamS("AbsOutput");

    if (s1[0] == 'T' || s1[0] == 't') {
      EFT_terms[ID]->abs_output = true;
      if ( EFT_terms[ID]->input_type == "fixed" )
	hf_errlog(23052302, "F: for fixed input, absolute output is not supported");
    }
  }

  if ( td->hasParam("NoCentral") ) {
    string s1 = td->getParamS("NoCentral");

    if (s1[0] == 'T' || s1[0] == 't')
      EFT_terms[ID]->no_central = true;
    else
      EFT_terms[ID]->no_central = false;
  }

  //------------------------------------------------------------------
  // read the filenames of coefficients
  EFT_terms[ID]->read_input();
}

//______________________________________________________________________________
// Main function to compute results at an iteration
void ReactionEFT::compute(TermData* td, valarray<double> &val, map<string, valarray<double> > &err)
{
  td->actualizeWrappers();

  EFTReader* EFT_term = EFT_terms[td->id];

  // X.S.: getParamD returns a pointer, is this address fixed?
  //vector<double> val_EFT_param;  
  valarray<double> val_EFT_param(100);  
  // val_EFT_param.reserve(EFT_term->num_param);

  int i=0;
  for (string EFT_param : EFT_term->name_EFT_param) {
    // val_EFT_param.push_back(*td->getParamD(EFT_param));
    val_EFT_param[i] = *td->getParamD(EFT_param);
    i++;
  }

  EFT_term->initIter(val_EFT_param);

  EFT_term->calcXSec(val);

  if (debug > 0) {
    std::cout << "=======================================================" << std::endl;
    std::cout << "ReactionEFT.compute: cross section: " << std::endl;
    for (double v : val){
      std::cout << v << ", ";
    }
    std::cout << std::endl;
  }
      
  // for ( std::size_t i=0; i<val.size(); i++) val[i] = (i<cs.size()) ? cs[i] : 0.;

}
