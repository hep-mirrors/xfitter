// DB 08/2017 
   /*
     @file ReactionEFT.cc
     @date 2016-12-06
     @author  AddReaction.py
     Created by  AddReaction.py on 2016-12-06
     Adapted from ReactionCIJET
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
  if ( td->hasParam("debug") )
    debug = td->getParamI("debug");

  // read the list of EFT parameters
  string list_EFT_param = td->getParamS("listEFTParam");
  // getNameEFTParam(list_EFT_param);

  if (list_EFT_param.size() == 0) {
    hf_errlog(23040601,"I: list of EFT parameters is empty");
  } else {
    stringstream ss(list_EFT_param);
    string param_name;
    while(getline(ss, param_name, ',')) {
      name_EFT_param.push_back(param_name);
    }
  }

  //------------------------------------------------------------------
  // read the filenames of coefficients
  if ( td->hasParam("Filename") ) {
    // --- Read file and instantiate EFT
    const std::string filename = td->getParamS("Filename");

    hf_errlog(23032801,"I: Reading EFT file "+filename);
    fname_list.push_back(filename);
    // EFT_terms.insert(std::make_pair(ID,vector<EFTReaction*>{new EFTReaction(filename, this)} ));
    EFT_terms.insert(std::make_pair(ID, new EFTReader(fname_list, debug, this) ));
  }
  else if ( td->hasParam("Filenames") ) {
    const std::string filenames = td->getParamS("Filenames");
    std::stringstream ss(filenames);
    std::string filename;
    while (std::getline(ss, filename, ',')) {
      hf_errlog(23032801,"I: Reading EFT file "+filename);
      fname_list.push_back(filename);
      // EFT_terms[ID].push_back(new EFTReaction(filename,this));
    }
    EFT_terms.insert(std::make_pair(ID, new EFTReader(fname_list, debug, this) ));
  }
  else {
    hf_errlog(23032803,"F:No EFT file specified. Please provide 'Filename' or 'Filenames' to reaction EFT.");
    return;
  }

  //------------------------------------------------------------------
  // read the coefficients
  EFT_terms[ID]->init(name_EFT_param); // read the coeff. from the input file
}

//______________________________________________________________________________
// Main function to compute results at an iteration
void ReactionEFT::compute(TermData* td, valarray<double> &val, map<string, valarray<double> > &err)
{
  td->actualizeWrappers();
  // getParamD returns a pointer, is this address fixed?
  vector<double> val_EFT_param;  
  val_EFT_param.reserve(name_EFT_param.size());

  for (string EFT_param : name_EFT_param) {
    val_EFT_param.push_back(*td->getParamD(EFT_param));
  }

  // vector<double> xsec_one_dataset;
  // xsec_one_dataset.reserve(val.size());
  int ID = td->id;

  EFTReader* EFT_term = EFT_terms[ID];
  EFT_term->setValEFT(val_EFT_param);
  std::vector<double> cs = EFT_term->calcxsec();

  if (debug > 0) {
    std::cout << "=======================================================" << std::endl;
    std::cout << "ReactionEFT.compute: cross section: " << std::endl;
    for (double v : cs){
      std::cout << v << ", ";
    }
    std::cout << std::endl;
  }
      
 
  if ( val.size() != cs.size() ){
    std::cout << "=======================================================" << std::endl;
    std::cout << val.size() << ", " <<  cs.size() << std::endl;
    std::cout << "=======================================================" << std::endl;
    hf_errlog(23032804,"F: Size of cross section array does not match data.");
  }
  for ( std::size_t i=0; i<val.size(); i++) val[i] = (i<cs.size()) ? cs[i] : 0.;

}
