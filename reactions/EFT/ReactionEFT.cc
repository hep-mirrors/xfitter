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
 
  int ID = td->id;  // ID=dataSetID -> updated to termdata td
  // read the list of EFT parameters
  string list_EFT_param = td->getParamS("listEFTParam");
  // getNameEFTParam(list_EFT_param);
  // a7: remove `[` and `]`
  list_EFT_param.erase(0, 1);
  list_EFT_param.erase(list_EFT_param.size()-1);

  stringstream ss(list_EFT_param);
  string temp;
  while(getline(ss, temp, ',')) { // todo: what if there is additional whitespace?
    name_EFT_param.push_back(temp);
  }

  //------------------------------------------------------------------
  // read the filenames of coefficients
  if ( td->hasParam("Filename") ) {
    // --- Read file and instantiate EFT
    const std::string filename = td->getParamS("Filename");

    hf_errlog(23032801,"I: Reading EFT file "+filename);
    EFT_all_dataset.insert(std::make_pair(ID,vector<EFTReaction*>{new EFTReaction(filename,this)} ));
  }
  else if ( td->hasParam("Filenames") ) {
    const std::string filenames = td->getParamS("Filenames");
    std::stringstream ss(filenames);
    std::string filename;
    while (std::getline(ss, filename, ',')) {
      hf_errlog(23032801,"I: Reading EFT file "+filename);
      EFT_all_dataset[ID].push_back(new EFTReaction(filename,this));
    }
  }
  else {
    hf_errlog(23032803,"F:No EFT file specified. Please provide 'Filename' or 'Filenames' to reaction EFT.");
    return;
  }

  //------------------------------------------------------------------
  // read the coefficients
  for ( EFTReaction* EFT_reader_one_file : EFT_all_dataset[ID] ) {
    // int num_bin = *td->getNBins(); //? TODO Does this really work? generally number of bins is file-dependent
    int num_bin = 7; // TODO !!
    EFT_reader_one_file->setinit(num_bin, name_EFT_param); //a7: read the coeff. from the input file
  }
}

//______________________________________________________________________________
// Main function to compute results at an iteration
void ReactionEFT::compute(TermData* td, valarray<double> &val, map<string, valarray<double> > &err)
{
  td->actualizeWrappers();
  // Get relevant parameters:
  // a7: read values of all EFT parameters; shall we try to save the time of memory allocation?
  // getParamD returns a pointer, is this address fixed?
  vector<double> val_EFT_param;  
  for (string EFT_param : name_EFT_param) {
    val_EFT_param.push_back(*td->getParamD(EFT_param));
  }

  vector<double> xsec_one_dataset;
  xsec_one_dataset.reserve(val.size());
  int ID = td->id;

  for ( EFTReaction* EFT_reader_one_file : EFT_all_dataset[ID] ) {
    EFT_reader_one_file->setValEFT(val_EFT_param);
    std::vector<double> cs = EFT_reader_one_file->calcxsec();
    // for ( std::size_t i=0; i<cs.size(); i++) xsec_one_dataset.push_back(cs[i]);
    for (double x: cs) xsec_one_dataset.push_back(x);
  }
 
  if ( val.size()!=xsec_one_dataset.size() )
    hf_errlog(23032804,"F: Size of EFT cross section array does not match data.");
  for ( std::size_t i=0; i<val.size(); i++) val[i] = (i<xsec_one_dataset.size()) ? xsec_one_dataset[i] : 0.;
}

//---------------------------------
// void getNameEFTParam(string list_EFT_param){
//   string s=list_EFT_param;
  
//   //Remove opening and closing brackets
//   s.erase(0, 1);
//   s.erase(input.size()-1);

//   stringstream ss(s);
//   string temp;
//   while(getline(ss, temp, ",")) {
//     name_EFT_Param.push_back(temp);
//   }
// };
