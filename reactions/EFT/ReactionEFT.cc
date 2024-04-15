   /*
     @file ReactionEFT.cc
     @date 2023-03
     @author X.M. Shen <xmshen137@gmail.com>
   */

#include "ReactionEFT.h"
#include <chrono>

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
  // int find_input_Q = 0;

  //------------------------------------------------------------------
  // construct EFTTerms
  /* 
  // allow for concatenation of bins in several files
  // 2024-04-10 X.M. Shen: practically never beneficial
  if ( td->hasParam("Filenames") ) {
    if (find_input_Q > 0) 
      hf_errlog(23052301, "F: there can be only one input in Filename/FileName/Filenames");
    find_input_Q = 1;

    const std::string filenames = td->getParamS("Filenames");
    std::stringstream ss(filenames);
    std::string filename;
    while (std::getline(ss, filename, ',')) {
      hf_errlog(23032801,"I: Reading EFT file "+filename);
      fname_list.push_back(filename);
    }
  }
  */
  std::cout << "==================================================================" << std::endl;
  std::cout << "EFT reaction: initialization for term " << ID << std::endl;

  if ( td->hasParam("FixedInput") ) {
    // find_input_Q = 1;
    const std::string filename = td->getParamS("FixedInput");
    hf_errlog(23032801,"I: Reading EFT file "+filename);
    fname_list.push_back(filename);
    EFT_terms.insert(std::make_pair(ID, new EFTTerm(fname_list, "fixed", this) ));

    std::cout << "EFT reaction: FixedInput found. AbsOutput is not supported in this mode." << std::endl;
  }
  else if ( td->hasParam("FileName") ) {
    // if (find_input_Q > 0) 
    //   hf_errlog(23052304, "F: expect only one of FileName");
    // find_input_Q = 2; // =1/2 for fixed/mixed input
    const std::string filename = td->getParamS("FileName");
    hf_errlog(23032801,"I: Reading EFT file "+filename);
    fname_list.push_back(filename);
    EFT_terms.insert(std::make_pair(ID, new EFTTerm(fname_list, "mixed", this) ));

    std::cout << "EFT reaction: FileName found" << std::endl;
  }
  else {
    hf_errlog(23032803,"F: EFT file should be specified with FileName in TermInfo");
  }

  //------------------------------------------------------------------
  // read the list of EFT parameter names
  if ( td->hasParam("ListEFTParam") ) {
    string list_EFT_param = td->getParamS("ListEFTParam");
    vector<string> name_EFT_param;

    if (list_EFT_param.size() == 0) {
      hf_errlog(23040601,"I: list of EFT parameters is empty. Not tested!");
    } else {
      stringstream ss(list_EFT_param);
      string param_name;
      while(getline(ss, param_name, ',')) {
	name_EFT_param.push_back(param_name);
      }
    }

    EFT_terms[ID]->initParamName(name_EFT_param);
    std::cout << "EFT reaction: ListEFTParam found" << std::endl;
  }
  else {
    hf_errlog(23061603,"F: ListEFTParam should be defined in TermInfo");
  }

  //-------------------------------------------------------
  // optional parameters; default values given in header file
  if ( td->hasParam("Debug") ) {
    EFT_terms[ID]->debug = td->getParamI("Debug");
    std::cout << "EFT reaction: set Debug" << std::endl;
  }

  if ( td->hasParam("AbsOutput") ) {
    string s1 = td->getParamS("AbsOutput");

    if (s1[0] == 'T' || s1[0] == 't') {
      EFT_terms[ID]->abs_output = true;
      if ( EFT_terms[ID]->input_type == "fixed" )
	hf_errlog(23052302, "F: for fixed input, absolute output is not supported");

      std::cout << "EFT reaction: set AbsOutput = true "  << std::endl;
      if ( td->hasParam("FixedInput") ) {
	hf_errlog(24041107,"S: EFT reaction: AbsOutput is not supported in FixedInput mode");
      }
    } 
    else
      std::cout << "EFT reaction: set AbsOutput = false " << std::endl;    
  }
  else
    std::cout << "EFT reaction: AbsOutput not set, = false by default" << std::endl;    

  if ( td->hasParam("NoCentral") ) {
    string s1 = td->getParamS("NoCentral");

    if (s1[0] == 'T' || s1[0] == 't') {
      EFT_terms[ID]->no_central = true;
      std::cout << "EFT reaction: set NoCentral = true" << std::endl; 
    }
    else {
      EFT_terms[ID]->no_central = false;
      std::cout << "EFT reaction: set NoCentral = false" << std::endl; 
    }
  }
  else
    std::cout << "EFT reaction: NoCentral not set, = false by default" << std::endl;    

  if ( td->hasParam("Norm") ) {
    string s1 = td->getParamS("Norm");
    if (s1[0] == 'T' || s1[0] == 't') {
      EFT_terms[ID]->normQ = true;
      std::cout << "EFT reaction: set Normalization = true" << std::endl; 
    }
    else
      std::cout << "EFT reaction: set Normalization = false" << std::endl; 
  }
  else
    std::cout << "EFT reaction: Norm not set, = false by default" << std::endl;    

  if ( td->hasParam("xiR") ) {
    EFT_terms[ID]->xi_ren = *(td->getParamD("xiR"));
    std::cout << "EFT reaction: set renormalization scale factor = " << EFT_terms[ID]->xi_ren << std::endl; 
  }
  else
    std::cout << "EFT reaction: renormalization scale factor not set, = 1.0 by default" << std::endl; 

  if ( td->hasParam("xiF") ) {
    EFT_terms[ID]->xi_fac = *(td->getParamD("xiF"));
    std::cout << "EFT reaction: set factorization scale factor = " << EFT_terms[ID]->xi_fac << std::endl; 
  }
  else
    std::cout << "EFT reaction: factorization scale factor not set, = 1.0 by default" << std::endl; 
  

  //------------------------------------------------------------------
  EFT_terms[ID]->readInput();

  std::cout << "EFT reaction: initialization completed for term " << ID << std::endl;
  std::cout << "======================================================="<< std::endl;
}



//______________________________________________________________________________
void ReactionEFT::freeTerm(TermData* td) {
  // todo
  // pineappl_grid_delete(data->grids[i]);
}



//______________________________________________________________________________
// Main function to compute results at an iteration
void ReactionEFT::compute(TermData* td, valarray<double> &val, map<string, valarray<double> > &err)
{
  auto startTime = std::chrono::high_resolution_clock::now();

  td->actualizeWrappers();

  EFTTerm* EFT_term = EFT_terms[td->id];
  valarray<double> val_EFT_param(100);  

  int i=0;
  for (string EFT_param : EFT_term->name_EFT_param) {
    val_EFT_param[i] = *td->getParamD(EFT_param);

    // if (EFT_term->debug > 2) { // already in initIter
    //   std::cout << "EFT_param: " << EFT_param << " = " << val_EFT_param[i] << std::endl;
    // }
    i++;
  }

  EFT_term->initIter(val_EFT_param);

  auto endTime1 = std::chrono::high_resolution_clock::now();

  EFT_term->calcXSec(val);

  if (EFT_term->debug > 0) {
    auto endTime2 = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(endTime1 - startTime).count(); // time cost for initialization: slow
    auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(endTime2 - endTime1).count();  // time cost for calc.: fast
    num_comp++; // number of computations = number of iteration * number of EFT terms in all data sets
    time_comp += duration1 + duration2;

    if (EFT_term->debug > 1) {
      std::cout << "=======================================================" << std::endl;
      std::cout << "ReactionEFT.compute[term=" << td->id << "]:" ;
      for (double v : val) std::cout << v << ", ";
      std::cout << std::endl;
    }

    std::cout << "Time cost(milliseconds)[term=" << td->id << "] [init, calc, total/num_comp]=" 
	      << duration1 << ", " 
	      << duration2 << ", " 
	      << time_comp / num_comp
	      << std::endl;
    std::cout << "=======================================================" << std::endl;
  } // end debug
      
}
