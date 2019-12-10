#include "xfitter_steer.h"
#include "xfitter_pars.h"
#include "xfitter_cpp_base.h"
#include"xfitter_cpp.h"

#include "BaseEvolution.h"
#include "BasePdfDecomposition.h"
#include "BasePdfParam.h"
#include "BaseMinimizer.h"
#include"ReactionTheory.h"
#include <dlfcn.h>
#include <iostream>
#include<fstream>
#include <yaml-cpp/yaml.h>
#include <Profiler.h>
using std::string;
using std::cerr;

extern std::map<string,string> gReactionLibs;

/*
\brief instantiate an object from a dynamically loaded library
\param moduleType
  moduleType is:
  "pdfparam"  for PDF parameterisations
  "pdfdecomp" for PDF decompositions
  "evolution" for evolutions
  "minimizer" for minimizers
  "reaction"  for reactions
\details
  Used to create evolutions, decompositions, parameterisations, reactions could be used to create minimizer
  
  The object is loaded using a (void*)create(instanceName) function loaded from dlopen-ed from an .so (shared object) library (module)
  The name of the loaded module library is
    "lib"+moduleType+className+".so"
  
  For moduleType in ["pdfparam", "pdfdecomp", "evolution"]
    instanceName is passed to create
  For moduleType in ["minimizer", "reaction"]
    create is called without arguments, and instanceName must be empty.
    (minimizers and reactions are supposed to be singletons)
*/
void* createDynamicObject(const string& moduleType, const string& className,const string& instanceName=""){

  using std::cerr;
  using std::endl;
  //On first call, initialize module prefix
  static bool first_call = true;
  //The following macro is defined in src/CMakeLists.txt
  //By default, it is INSTALL_PREFIX/lib/xfitter/
  static string module_prefix = XFITTER_DEFAULT_MODULE_PATH;
  const char* XFITTER_MODULE_PATH_ENV_NAME="XFITTER_MODULE_PATH";
  if (first_call) {
    first_call = false;
    const char* from_env = getenv(XFITTER_MODULE_PATH_ENV_NAME);
    if (from_env != nullptr) {
      module_prefix = from_env;
    }

    //Make fure module_prefix ends in "/"
    if ( module_prefix.empty() ) {
      hf_errlog(19081600,"W: XFITTER_MODULE_PREFIX environment variable is defined and empty, will search for dynamically loaded libraries in the working directory");
    } else if( module_prefix.back() != '/' ) {
      module_prefix += '/';
    }
  }

  //form the path to the loaded library
  string libpath = module_prefix + "lib" + moduleType + className + ".so";
  //load the library
  void* shared_library = dlopen(libpath.c_str(), RTLD_NOW);
  //by the way, do we ever call dlclose? I don't think so... Maybe we should call it eventually. --Ivan Novikov

  //error if failed to load library
  if ( shared_library == nullptr ){
    cerr<<"[ERROR] dlopen() error while trying to open shared library for class \""<<className<<"\":\n"
      <<dlerror()<<"\n"
      "xFitter failed to load module "<<libpath<<
      "\nMake sure that the class name \""<<className<<"\" in the YAML steering is correct, and that the required module is installed\n"
      "If your modules directory is located elsewhere, set the environment variable "<<XFITTER_MODULE_PATH_ENV_NAME<<" to the correct directory"
      "\n[/ERROR]"<<endl;
    hf_errlog(18091900,"F: dlopen() error, see stderr");
  }
  //reset errors
  dlerror();

  //load the create() function from the library
  void* create = dlsym(shared_library,"create");
  if ( create == nullptr ){
    cerr<<"[ERROR] dlsym() failed to find \"create\" function for class \""<<className<<"\":\n"<<dlerror()<<"in loaded module "<<libpath<<"\n[/ERROR]"<<endl;
    hf_errlog(18091902,"F:dlsym() error, see stderr");
  }

  void* obj;
  if( moduleType=="reaction" or moduleType=="minimizer"){
    if (not instanceName.empty()) {
      cerr<<"[ERROR] "<<__func__<<"() called with invalid arguments:\n"
      "moduleType=\""<<moduleType<<"\n"
      "className=\""<<className<<"\" (should have been empty)\n"
      "instanceName=\""<<instanceName<<"\" (should have been empty =\"\")\n"
      "a "<<moduleType<<" should be a singleton, it cannot have an instanceName."
      "Somebody go fix the code"<<endl;
      hf_errlog(19081601,"F: Tried to name a create a named singleton, see std");
    }
    //call create without arguments
    obj=((void*(*)())create)();
  }else{
    //pass instanceName to create
    obj=((void*(*)(const char*))create)(instanceName.c_str());
  }

  //Name consistency check: the requested name and the one reported by the class must match
  string reportedName;

  //get reported name depending on base class
  if( moduleType == "pdfdecomp" ){
    reportedName=((xfitter::BasePdfDecomposition*)obj)->getClassName();
  }else if( moduleType == "evolution" ){
    reportedName=((xfitter::BaseEvolution*)obj)->getClassName();
  }else if( moduleType == "reaction" ){
    reportedName=((ReactionTheory*)obj)->getReactionName();
  }else{
    //no check for pdfparam or minimizer, just return
    return obj;
  }
  if( reportedName != className ) {
    cerr<<"[ERROR] class name mismatch:\n"
      "\""<<className<<" expected\n"
      "\""<<reportedName<<" reported by the class\n"
      "for dynamically loaded module "<<libpath<<endl;
    hf_errlog(19081602,"F: Class name check failed, see stderr");
  }
  return obj;
}

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
    hf_errlog(18082950,s.str().c_str());
  }
  string classname=classnameNode.as<string>();
  BaseEvolution* evolution = (BaseEvolution*)createDynamicObject("evolution", classname, name);
  //Note that unlike in the pervious version of this function, we do not set decompositions for evolutions
  //Evolution objects are expected to get their decomposition themselves based on YAML parameters, during atStart
  try{
    evolution->atStart();
  }catch(const YAML::Exception&){
    cerr<<"[ERROR] Unhandled YAML exception during initialization of evolution \""<<name<<"\". Check that the provided parameters are correct. Rethrowing the exception..."<<endl;
    throw;
  }
  // Store the newly created evolution on the global map
  XFITTER_PARS::gEvolutions[name] = evolution;
  return evolution;
}
BasePdfDecomposition*get_pdfDecomposition(string name){
  try{
    if(name=="")name=XFITTER_PARS::getDefaultDecompositionName();
    auto it=XFITTER_PARS::gPdfDecompositions.find(name);
    if(it!=XFITTER_PARS::gPdfDecompositions.end())return it->second;
    string classname = XFITTER_PARS::getDecompositionNode(name)["class"].as<string>();
    BasePdfDecomposition* ret = (BasePdfDecomposition*)createDynamicObject("pdfdecomp", classname, name);
    ret->atStart();
    XFITTER_PARS::gPdfDecompositions[name]=ret;
    return ret;
  }catch(const YAML::TypedBadConversion<string>&ex){
    using namespace std;
    YAML::Node node = XFITTER_PARS::getDecompositionNode(name)["class"];
    if (!node){
      cerr<<"[ERROR] No class specified for decomposition \""<<name<<"\": missing node Decompositions/"<<name<<"/class"<<endl;
      hf_errlog(18092401, "F: No class specified for a decomposition, see stderr");
    }
    cerr<<"[ERROR] Failed to convert class specification for decomposition \""<<name<<"\" to string, node Decompositions/"<<name<<"/class"<<endl;
    cerr<<"Trying to print the node:"<<endl;
    cerr<<node<<endl;
    cerr<<"[/ERROR]"<<endl;
    hf_errlog(18092402, "F: Bad class specification for a decomposition, see stderr");
  }
}
BasePdfParam*getParameterisation(const string&name){
  try{
    auto it=XFITTER_PARS::gParameterisations.find(name);
    if(it!=XFITTER_PARS::gParameterisations.end())return it->second;
    //Else create a new instance
    string classname=XFITTER_PARS::getParameterisationNode(name)["class"].as<string>();
    BasePdfParam* ret = (BasePdfParam*)createDynamicObject("pdfparam", classname, name);
    ret->atStart();
    XFITTER_PARS::gParameterisations[name]=ret;
    return ret;
  }catch(YAML::InvalidNode&ex){
    const int errcode=18092400;
    const char*errmsg="F: YAML::InvalidNode exception while creating parameterisation, details written to stderr";
    using namespace std;
    cerr<<"[ERROR]"<<__func__<<'('<<name<<')'<<endl;
    YAML::Node node=XFITTER_PARS::getParameterisationNode(name);
    if(!node.IsMap()){
      cerr<<"Invalid node Parameterisations/"<<name<<"\nnode is not a map\n[/ERROR]"<<endl;
      hf_errlog(errcode,errmsg);
    }
    YAML::Node node_class=node["class"];
    if(!node_class.IsScalar()){
      if(node_class.IsNull())cerr<<"Missing node Parameterisations/"<<name<<"/class";
      else cerr<<"Invalid node Parameterisations/"<<name<<"/class\nnode is not a scalar";
      cerr<<"\n[/ERROR]"<<endl;
      hf_errlog(errcode,errmsg);
    }
    cerr<<"Unexpected YAML exception\nNode:\n"<<node<<"\n[/ERROR]"<<endl;
    hf_errlog(errcode,errmsg);
  }
}


BaseMinimizer* get_minimizer() {
  std::string name = XFITTER_PARS::getParamS("Minimizer");

  // Check if already present
  if (XFITTER_PARS::gMinimizer != nullptr ) {
    return  XFITTER_PARS::gMinimizer;  //already loaded
  }

  // else load, initialize and return
  XFITTER_PARS::gMinimizer =(BaseMinimizer*) createDynamicObject("minimizer", name);
  XFITTER_PARS::gMinimizer->atStart();
  return XFITTER_PARS::gMinimizer;
}

ReactionTheory* getReaction(const string& name){
  //if already exists, return it
  auto it = gNameReaction.find(name);
  if ( it != gNameReaction.end() ) return it->second;
  //else create and return
  ReactionTheory* rt=(ReactionTheory*)createDynamicObject("reaction", name);
  gNameReaction[name] = rt;
  //initialize
  rt->atStart();
  return rt;
}

}


/// Temporary interface for fortran
extern "C" {
void init_minimizer_();
void run_minimizer_();
void report_convergence_status_();
void run_error_analysis_();
}

namespace xfitter{
  BaseEvolution*defaultEvolution=nullptr;//declared in xfitter_steer.h

  BaseEvolution* defaultEvolutionInstance() {
    if ( defaultEvolution==nullptr) defaultEvolution=xfitter::get_evolution();
    return defaultEvolution;
  }
}

void init_minimizer_() {
  /// atStart is called inside
  auto mini = xfitter::get_minimizer();
}

void run_minimizer_() {
  auto mini = xfitter::get_minimizer();
  /// get profiler too
  auto *prof = new xfitter::Profiler();

  prof->doProfiling();

  mini->doMinimization();
}

void report_convergence_status_(){
  //Get a status code from current minimizer and log a message, write status to file Status.out
  using namespace xfitter;
  auto status=get_minimizer()->convergenceStatus();
  //Write status to Status.out
  {
    std::ofstream f;
    f.open(stringFromFortran(coutdirname_.outdirname,sizeof(coutdirname_.outdirname))+"/Status.out");
    if(!f.is_open()){
      hf_errlog(16042807,"W: Failed to open Status.out for writing");
      return;
    }
    if(status==ConvergenceStatus::SUCCESS)f<<"OK";
    else f<<"Failed";
    f.close();
  }
  //Log status message
  switch(status){
    case ConvergenceStatus::NORUN:
      hf_errlog(16042801,"I: No minimization has run");
      break;
    case ConvergenceStatus::INACCURATE:
      hf_errlog(16042803,"S: Error matrix not accurate");
      break;
    case ConvergenceStatus::FORCED_POSITIVE:
      hf_errlog(16042804,"S: Error matrix forced positive");
      break;
    case ConvergenceStatus::SUCCESS:
      hf_errlog(16042802,"I: Fit converged");
      break;
    case ConvergenceStatus::NO_CONVERGENCE:
      hf_errlog(16042805,"S: No convergence");
      break;
    case ConvergenceStatus::ERROR:
      hf_errlog(16042806,"F: Minimizer error");
      break;
  }
}

void run_error_analysis_() {
  auto mini = xfitter::get_minimizer();
  mini->errorAnalysis();
}

namespace xfitter{
void updateAtConfigurationChange(){
  //Call atConfigurationChange for each evolution and for each decomposition
  for(map<string,BaseEvolution*>::const_iterator it=XFITTER_PARS::gEvolutions.begin();it!=XFITTER_PARS::gEvolutions.end();++it){
    it->second->atConfigurationChange();
  }
  for(map<string,BasePdfDecomposition*>::const_iterator it=XFITTER_PARS::gPdfDecompositions.begin();it!=XFITTER_PARS::gPdfDecompositions.end();++it){
    it->second->atConfigurationChange();
  }
}
}
