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
using std::string;

extern std::map<string,string> gReactionLibs;

namespace xfitter {

  BaseEvolution*get_evolution(string name){
    if(name=="")name=XFITTER_PARS::getDefaultEvolutionName();
    // Check if already present
    if(XFITTER_PARS::gEvolutions.count(name)==1){
      return XFITTER_PARS::gEvolutions.at(name);
    }
    // Else create a new instance of evolution
    YAML::Node instanceNode=XFITTER_PARS::getEvolutionNode(name);
    YAML::Node classnameNode=instanceNode["class"];
    if(!classnameNode.IsScalar()){
      std::ostringstream s;
      s<<"F:Failed to get evolution \""<<name<<"\": evolution must have a node \"class\" with the class name as a string";
      hf_errlog(18082902,s.str().c_str());
    }
    string classname=classnameNode.as<string>();
    string libname;
    {auto it=gReactionLibs.find(classname);
    if(it==gReactionLibs.end()){
      std::ostringstream s;
      s<<"F: Failed to get evolution \""<<name<<"\": unknown evolution class name \""<<classname<<"\"";
      hf_errlog(18082903,s.str().c_str());
    }
    libname=it->second;
    }
    // load the library:
    //dlopen may be called multiple times for the same library and should return the same handle each time
    void*shared_library=dlopen((PREFIX+string("/lib/")+libname).c_str(),RTLD_NOW);
    //By the way, do we ever call dlclose? I don't think so... Maybe we should call it eventually. --Ivan Novikov
    if(shared_library==NULL){ 
      std::cout<<dlerror()<<std::endl;      
      hf_errlog(18071303,"F: Evolution shared library ./lib/"  + libname  +  " not present for evolution" + name + ". Check Reactions.txt file");
    }
    // reset errors
    dlerror();

    BaseEvolution*evolution=((create_evolution*)dlsym(shared_library,"create"))(name.c_str());

    //Note that unlike in the pervious version of this function, we do not set decompositions for evolutions
    //Evolution objects are expected to get their decomposition themselves based on YAML parameters, during initFromYaml

    evolution->initFromYaml(instanceNode);
    // Store the newly created evolution on the global map
    XFITTER_PARS::gEvolutions[name] = evolution;
    return evolution;
  }
  BasePdfDecomposition*get_pdfDecomposition(string name){
    if(name=="")name=XFITTER_PARS::getDefaultDecompositionName();
    auto it=XFITTER_PARS::gPdfDecompositions.find(name);
    if(it!=XFITTER_PARS::gPdfDecompositions.end())return it->second;
    // Load corresponding shared library:
    //Note: apparently we store all shared lib names in the map gReactionLibs, whether the libs are for reaction or whatever else. This naming is misleading.
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
  //TODO: reimplement for new interface with multiple evolutions
  //auto evol = xfitter::get_evolution();
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

