#include "xfitter_steer.h"
#include "xfitter_pars.h"
#include "xfitter_cpp_base.h"

#include "BaseEvolution.h"
#include <dlfcn.h>
#include <iostream>

extern std::map<string,string> gReactionLibs;

namespace xfitter {

  void load_evolution(std::string name) {
    if (name == "") {
      // get the name from the map
      name = XFITTER_PARS::getParameterS("Evolution");
    }    
    // Load corresponding shared library:
    string libname = gReactionLibs[name];
    if ( libname == "") {
      hf_errlog(18071302,"F: Shared library for evolution "+name+" not found");
    }

    // load the library:
    void *evolution_handler = dlopen((PREFIX+string("/lib/")+libname).c_str(), RTLD_NOW);
    if (evolution_handler == NULL)  { 
      std::cout  << dlerror() << std::endl;      
      hf_errlog(18071303,"F: Evolution shared library ./lib/"  + libname  +  " not present for evolution" + name + ". Check Reactions.txt file");
    }

     // reset errors
    dlerror();

    create_evolution *dispatch_ev = (create_evolution*) dlsym(evolution_handler, "create");
   
    BaseEvolution *evolution = dispatch_ev();

    /// now we want to store it somewhere XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
  }
  void load_pdfDecomposition(std::string name) {
    if (name == "") {
      // get the name from the map
      name = XFITTER_PARS::getParameterS("PDFDecomposition");
    }    
    // Load corresponding shared library:
    string libname = gReactionLibs[name];
    if ( libname == "") {
      hf_errlog(19072302,"F: Shared library for pdf decomposition "+name+" not found");
    }

    // load the library:
    void *pdfDecomposition_handler = dlopen((PREFIX+string("/lib/")+libname).c_str(), RTLD_NOW);
    if (pdfDecomposition_handler == NULL)  { 
      std::cout  << dlerror() << std::endl;      
      hf_errlog(18071303,"F: PdfDecomposition shared library ./lib/"  + libname  +  " not present for pdfDecomposition" + name + ". Check Reactions.txt file");
    }
  
  }
}
