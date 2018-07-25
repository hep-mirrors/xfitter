#include "xfitter_steer.h"
#include "xfitter_pars.h"
#include "xfitter_cpp_base.h"

#include "BaseEvolution.h"
#include "BasePdfDecomposition.h"
#include <dlfcn.h>
#include <iostream>
#include <yaml-cpp/yaml.h>

extern std::map<string,string> gReactionLibs;

namespace xfitter {

  BaseEvolution* get_evolution(std::string name) {
    if (name == "") {
      // get the name from the map
      name = XFITTER_PARS::getParameterS("Evolution");
    }
    
    // Check if already present
    if (XFITTER_PARS::gEvolutions.count(name) == 1) {
      return XFITTER_PARS::gEvolutions[name];  //already loaded
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

    // Now we attach corresponding PDFdecomposition. First try specific, next: global
    std::string pdfDecomp =  XFITTER_PARS::getParameterS("PDFDecomposition");

    // XXXXXXXXXXXXXXXXXXXXXXXX
    if ( XFITTER_PARS::gParametersY.count(name) > 0 ) {
      auto evolNode = XFITTER_PARS::gParametersY[name];
      if ( evolNode["PDFDecomposition"] ) {
	pdfDecomp = evolNode["PDFDecomposition"].as<std::string>();
	std::cout << " here here \n";
      }
      else {
	std::cout << " ho here \n";	
      }
    }

    std::cout << "PDF decomp=" << pdfDecomp << "\n";

    // Get corresponding PDF decomposition:
    BasePdfDecomposition* pdfD = get_pdfDecomposition(pdfDecomp);
    evolution->SetPdfDecomposition( pdfD->f0() );
    // Store on the map
    XFITTER_PARS::gEvolutions[name] = evolution;
    return evolution;
  }

  BasePdfDecomposition* get_pdfDecomposition(std::string name) {
    if (name == "") {
      // get the name from the map
      name = XFITTER_PARS::getParameterS("PDFDecomposition");
    }
    // Check if already present
    if (XFITTER_PARS::gPdfDecompositions.count(name) == 1) {
      std::cout << " already present " << std::endl;
      return  XFITTER_PARS::gPdfDecompositions[name];  //already loaded
    }
    
    // Load corresponding shared library:
    string libname = gReactionLibs[name];
    if ( libname == "") {
      hf_errlog(18072302,"F: Shared library for pdf decomposition "+name+" not found");
    }

    // load the library:
    void *pdfDecomposition_handler = dlopen((PREFIX+string("/lib/")+libname).c_str(), RTLD_NOW);
    if (pdfDecomposition_handler == NULL)  { 
      std::cout  << dlerror() << std::endl;      
      hf_errlog(18072303,"F: PdfDecomposition shared library ./lib/"  + libname  +  " not present for pdfDecomposition" + name + ". Check Reactions.txt file");
    }

         // reset errors
    dlerror();

    create_pdfDecomposition *dispatch_decomp = (create_pdfDecomposition*) dlsym(pdfDecomposition_handler, "create");
    BasePdfDecomposition *pdfDecomp = dispatch_decomp();
    pdfDecomp->initAtStart("");

    // store on the map
    XFITTER_PARS::gPdfDecompositions[name] = pdfDecomp;

    return pdfDecomp;
  }
}
