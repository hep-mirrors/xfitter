#include "xfitter_steer.h"
#include "xfitter_pars.h"
#include "xfitter_cpp_base.h"

#include "BaseEvolution.h"
#include "BasePdfDecomposition.h"
#include "BaseMinimizer.h"
#include <dlfcn.h>
#include <iostream>
#include <yaml-cpp/yaml.h>
#include <Profiler.h>

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

    // Init it:
    evolution->initAtStart();
    
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


  BaseMinimizer* get_minimizer() {
    std::string name = XFITTER_PARS::getParameterS("Minimizer");
    
    // Check if already present
    if (XFITTER_PARS::gMinimizer != nullptr ) {
      return  XFITTER_PARS::gMinimizer;  //already loaded
    }
    
    // Load corresponding shared library:
    string libname = gReactionLibs[name];
    if ( libname == "") {
      hf_errlog(18081701,"F: Shared library for minimizer "+name+" not found");
    }

    // load the library:
    void *lib_handler = dlopen((PREFIX+string("/lib/")+libname).c_str(), RTLD_NOW);
    if ( lib_handler == nullptr )  { 
      std::cout  << dlerror() << std::endl;      
      hf_errlog(18081702,"F: Minimizer shared library ./lib/"  + libname  +  " not present for minimizer" + name + ". Check Reactions.txt file");
    }

         // reset errors
    dlerror();

    create_minimizer *dispatch_minimizer = (create_minimizer*) dlsym(lib_handler, "create");
    BaseMinimizer *minimizer = dispatch_minimizer();
    minimizer->initAtStart();

    // store on the map
    XFITTER_PARS::gMinimizer = minimizer;

    return XFITTER_PARS::gMinimizer;
  }

}


/// Temporary interface for fortran
extern "C" {
  void init_evolution_(); 
  void init_minimizer_();
  void run_minimizer_();
  void run_error_analysis_();
}

void init_evolution_() {
  auto evol = xfitter::get_evolution();
}

void init_minimizer_() {
  /// initAtStart is called inside
  auto mini = xfitter::get_minimizer();
}

void run_minimizer_() {
  auto mini = xfitter::get_minimizer();
  /// get profiler too
  auto *prof = new xfitter::Profiler();

  prof->doProfiling();
  
  mini->doMimimization();    
}

void run_error_analysis_() {
  auto mini = xfitter::get_minimizer();
  mini->errorAnalysis();    
}

