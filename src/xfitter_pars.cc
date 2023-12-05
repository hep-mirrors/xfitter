/*!
 @file xfitter_pars.cc
 @date Sun 16 April 2017
 @author SG

 Contains functions to read parameters.yaml,
 global maps to store parameters,  and fortran interface functions.
 */

#include "xfitter_pars.h"
#include "xfitter_steer.h"
#include "xfitter_cpp.h"
#include "xfitter_cpp_base.h"
#include <string.h>
#include <cmath>
#include "BaseEvolution.h"
#include "BasePdfParam.h"
#include "BasePdfDecomposition.h"
#include "BaseMinimizer.h"
#include"ReactionTheory.h"
#include <dlfcn.h>
#include "dependent_pars.h"
#include "ansi_codes.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

// Fortran bindings:
extern "C" {
  void parse_params_();  //!< Fortran callable parsing
  /// Interface to minuit parameters
  void addexternalparam_(const char name[],  const double &val,
       const double  &step,
       const double &min, const double &max,
       const double &prior, const double &priorUnc,
       const int &add,
       map<std::string,double*> *map,
       int len);
  void add_to_param_map_(map<std::string,double*>* map,double &value, int& global, char *name, int len);
  // Get parameters in fortran, for backward compatibility:
  double getparamd_(const char* name, int len);
  int getparami_(const char* name, int len);
  // Update of EWK/QCD parameters, can be fitted at each iteration.
  void update_pars_fortran_();
}

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

  // Check if the library with extension .so exists (this should be the case on Linux)
  struct stat info;
  int ret = stat(libpath.c_str(), &info);
  // If not try with extension .dylib (this should be the case on MacOS)
  if (ret != 0)
    libpath = module_prefix + "lib" + moduleType + className + ".dylib";

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

namespace XFITTER_PARS {

  // Global vars:
  xfitter::BaseMinimizer* gMinimizer(nullptr);
  YAML::Node rootNode;
  map <string, double*> gParameters;
  map <string, int>    gParametersI;
  map <string, string> gParametersS;
  map <string, vector<double> > gParametersV; ///< Vectors of double parameters
  map <string, vector<string> > gParametersVS; ///< Vectors of string (not double) parameters
  map <string, YAML::Node > gParametersY;      ///< Store complete nodes for complex cases

  // Also keep list of loaded evolutions here:
  map<string,xfitter::BaseEvolution*> gEvolutions;
  // Also keep list of loaded decompositions here:
  map<string,xfitter::BasePdfDecomposition*> gPdfDecompositions;
  // Also keep list of loaded parameterisations here:
  map<string,xfitter::BasePdfParam*> gParameterisations;

  // Functions to get parameters from corresponding maps but with better reporting of errors
  double*getParamD(const string&name){
    try{return gParameters.at(name);}
    catch(std::out_of_range&ex){
      cerr<<"[ERROR] Double parameter \""<<name<<"\" does not exist"<<flush;
      //Check for a potential incorrect definition of a fit parameter
      YAML::Node node = rootNode[name];
      if (node.IsSequence() or node.IsMap() and node["value"]) {
        cerr<<", but there is something that looks like a definition of this parameter at global scape."
        " It is possible that you wanted to fit this parameter, but did not define it correctly."
        "\nAll fitted parameters must be provided in YAML steering under the \"Parameters\" node and not in any other place."
        "\nBefore fitting theory parameters, it is also important to make sure that the used evolutions"
        " and reaction modules re-read the updated parameter value at each iteration, and not only at initialization."
        "\nxFitter will now crash";
      }
      cerr<<"; rethrowing out_of_range"<<endl;
      throw ex;//rethrow exception: makes it easier to debug who tried to get parameter
    }
  }
  int getParamI(const string&name){
    try{return gParametersI.at(name);}
    catch(std::out_of_range&ex){
      cerr<<"[ERROR] Int parameter \""<<name<<"\" does not exist; rethrowing out_of_range"<<endl;
      throw ex;//rethrow exception: makes it easier to debug who tried to get parameter
    }
  }
  string getParamS(const string&name){
    try{return gParametersS.at(name);}
    catch(std::out_of_range&ex){
      cerr<<"[ERROR] String parameter \""<<name<<"\" does not exist; rethrowing out_of_range"<<endl;
      throw ex;//rethrow exception: makes it easier to debug who tried to get parameter
    }
  }

  using namespace xfitter;
  xfitter::BasePdfDecomposition*getInputDecomposition(const YAML::Node&rootNode){
    YAML::Node node=rootNode["decomposition"];
    if(!node)return xfitter::get_pdfDecomposition("");//return default decomposition
    string name;
    try{
      name=node.as<string>();
    }catch(YAML::TypedBadConversion<string>ex){
      ostringstream s;s<<"W: YAML exception: "<<ex.what()<<"; while trying to extract decomposition name from node: "<<node<<"; using default decomposition name";
      hf_errlog(18082930,s.str());
      name=getDefaultDecompositionName();
    }
    if (name == "None")
      return nullptr;
    return xfitter::get_pdfDecomposition(name);
  }
  YAML::Node getEvolutionNode(const std::string&name){
    auto it=gParametersY.find("Evolutions");
    if(it==gParametersY.end()){
      hf_errlog(18082900,"F:Failed to get evolution "+name+": Evolutions node not found in parameters.yaml");
    }
    YAML::Node instanceNode=it->second[name];
    if(!instanceNode.IsMap()){
      std::ostringstream s;
      s<<"F:Failed to get evolution \""<<name<<"\": ";
      if(!instanceNode)s<<"no subnode with this name under the node Evolutions";
      else s<<"corresponding subnode is not of type Map";
      hf_errlog(18082901,s.str().c_str());
    }
    return instanceNode;
  }
  YAML::Node getDecompositionNode(const std::string&name){
    auto it=gParametersY.find("Decompositions");
    if(it==gParametersY.end()){
      hf_errlog(18082902,"F:Failed to get decomposition "+name+": \"Decompositions\" node not found");
    }
    YAML::Node instanceNode=it->second[name];
    if(!instanceNode.IsMap()){
      std::ostringstream s;
      s<<"F:Failed to get decomposition \""<<name<<"\": ";
      if(!instanceNode)s<<"no subnode with this name under the node Decompositions";
      else s<<"corresponding subnode is not of type Map";
      hf_errlog(18082903,s.str().c_str());
    }
    return instanceNode;
  }
  YAML::Node getParameterisationNode(const std::string&name){
    auto it=gParametersY.find("Parameterisations");
    if(it==gParametersY.end()){
      hf_errlog(18082904,"F:Failed to get parameterisation "+name+": \"Parameterisations\" node not found");
    }
    YAML::Node instanceNode=it->second[name];
    if(!instanceNode.IsMap()){
      std::ostringstream s;
      s<<"F:Failed to get parameterisation \""<<name<<"\": ";
      if(!instanceNode)s<<"no subnode with this name under the node Parameterisations";
      else s<<"corresponding subnode is not of type Map";
      hf_errlog(18082905,s.str().c_str());
    }
    return instanceNode;
  }

  string getDefaultEvolutionName(){
    auto it=gParametersS.find("DefaultEvolution");
    if(it!=gParametersS.end()) return it->second;
    return "default";
  }
  string getDefaultDecompositionName(){
    auto it=gParametersS.find("DefaultDecomposition");
    if(it!=gParametersS.end())return it->second;
    return "default";
  }

double getEvolutionParamD(const string& evName,const string& parName){
  YAML::Node parNode = getEvolutionNode(evName)[parName];
  if (!parNode){
    cerr<<"[ERROR] Missing parameter \""<<parName<<"\" for evolution \""<<evName<<"\""<<endl;
    hf_errlog(19051630, "F: Missing evolution parameter, see stderr");
  }
  try{
    return parNode.as<double>();
  }catch (YAML::TypedBadConversion<double>& ex){
    cerr<<"[ERROR] Failed to convert to double parameter \""<<parName<<"\" of evolution \""<<evName<<"\""<<endl;
    hf_errlog(19051631, "F: Failed to convert evolution parameter to double, see stderr");
  }
  std::abort();//unreachable
}

YAML::Node loadYamlFile(const string&filename){
  YAML::Node node;
  try{
    node=YAML::LoadFile(filename);
  }catch(const YAML::BadFile&ex){
    cerr<<"[ERROR] Failed to open yaml file "<<filename<<endl;
    hf_errlog(18092600,"F: Failed to open yaml file, see stderr");
  }
  return node;
}
/*
\brief compare two nodes, returns true if they are equal (works only if they are scalars or sequences, does not work for maps: returns false)
*/
bool areEqual(const YAML::Node&node1,const YAML::Node&node2,unsigned int recursionLimit=256){
    if(recursionLimit==0){
        hf_errlog(20032501,"F: Recursion limit reached while checking YAML list equality");
    }
    if(node1.IsScalar() && node2.IsScalar()){
        return (node1.as<string>() == node2.as<string>());
    }
    else if(node1.IsSequence() && node2.IsSequence()){
        for(size_t i = 0; i < node1.size(); i++){
            auto ret = areEqual(node1[i], node2[i], recursionLimit-1);
            if (!ret) return false;
        }
    }
    // otherwise we assume nodes are different: it is not clear how to compare maps
    return false;
}
/*
\brief Process "? !include" directives in the YAML steering
\details
  For each "? !include" directive in the loaded YAML tree,
  open the included file and insert its contents in place of the !include,
  possibly ignoring some keys if they are already defined locally.

  When looking for the included file, search relative to current working directory first.
  If the included file is not found there,
  AND the given path is not absolute (absolute paths begin with '/'),
  AND the given path is not explicitly relative to working dir (begin with '.'),
  then look in the system directory for standard xfitter YAML files,
  which is XFITTER_YAML_PATH=INSTALL_PREFIX/share/xfitter/

  This function operates on YAML maps, looking for "include" tag on key of each entry.

  Include expansion is recursive, which means that "!include" statements will be expanded in sub-maps too
  (on any indentation level of the YAML steering).

  Included files can also contain "!include" directives

  The recursionLimit is meant to protect from circular includes.
*/
void expandIncludes(YAML::Node&node,unsigned int recursionLimit=256){
  if(recursionLimit==0){
    hf_errlog(18092605,"F: Recursion limit reached while handling includes");
  }
  if(!node.IsMap())return;//maybe this should even be an error
  vector<YAML::Node>include_keys;
  for(YAML::iterator it=node.begin();it!=node.end();++it){
    YAML::Node key=it->first;
    YAML::Node val=it->second;
    if(key.Tag()=="!include"){
      if(!it->second.IsNull()){
        cerr<<"[ERROR] Value given after include key "<<key<<" (make sure there is no \":\" after the filename)"<<endl;
        hf_errlog(18092602,"F: Value after include key, see stderr");
      }
      include_keys.push_back(key);
    }else if(val.IsMap()){
      expandIncludes(val,recursionLimit-1);
    }
  }
  //Load and merge
  for(vector<YAML::Node>::const_iterator kit=include_keys.begin();kit!=include_keys.end();++kit){
    node.remove(*kit);
    string filename=(*kit).as<string>("");
    if(filename==""){
      cerr<<"[ERROR] Failed to parse include filename, node:\n"<<node<<"\n[/ERROR]"<<endl;
      hf_errlog(18092601,"F: Failed to parse include filename, see stderr");
    }
    //Search for include files first relative to current working directory, then in PREFIX/yaml
    if(!fileExists(filename)){
      //if filename starts with '/', it is an absolute path
      //if filename starts with '.', it is a path relative to current directory
      //In these two cases, do not search in PREFIX/yaml
      char c=filename[0];
      bool file_not_found=(c=='/'||c=='.');
      if(!file_not_found){
        string prefix_filename=XFITTER_YAML_PATH+filename;
        if(fileExists(prefix_filename))filename=prefix_filename;
        else file_not_found=true;
      }
      if(file_not_found){
        cerr<<"[ERROR] YAML include file "<<filename<<" not found\n"
        "Default search prefix is "<<XFITTER_YAML_PATH<<'\n'<<endl;
        hf_errlog(19040135,"F: YAML include file not found, see stderr");
      }
    }
    YAML::Node include=loadYamlFile(filename);
    if(!include.IsMap()){
      cerr<<"[ERROR] Root node in included file "<<filename<<" is not a map"<<endl;
      hf_errlog(18092603,"F: Trying to include something other than a map, see stderr");
    }
    expandIncludes(include,recursionLimit-1);
    for(YAML::const_iterator it=include.begin();it!=include.end();++it){
      if(it->first.IsScalar()){
        string key=it->first.Scalar();
        if(node[key]){
            if(!areEqual(it->second, node[key])){
              cout<<"[INFO] Option "<<key<<"="<<it->second<<" included from file "<<filename<<" is overridden by locally defined option "<<key<<"="<<node[key]<<endl;
              hf_errlog(18092604,"I: locally defined setting overrides included, see stdlog");
            }
        }else node[key]=it->second;
      }else{
        hf_errlog(19033101,"W: including YAML maps with non-scalar keys is poorly supported and does not handle overriding keys, use at your own risk!");
        node[it->first]=it->second;
      }
    }
  }
}
  // Parse @param node and return maps
  void parse_node(const YAML::Node& node,
                  std::map<string,double*>& dMap,
                  std::map<string,int>& iMap,
                  std::map<string,string>& sMap,
                  std::map<string,vector<double> >& vMap,
		  std::map<string,vector<string> >& vsMap,
                  std::map<string,YAML::Node> & yMap ){
    for ( YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
      YAML::Node key = it->first;
      YAML::Node value = it->second;
      // parameter name
      string p_name = key.as<string>();
      if (value.IsScalar()) {
        // Alright, store directly
        // Try to read as int, float, string:
        try {
          int i = value.as<int>();
          iMap[p_name] = i;
          continue;
        }
        catch (const std::exception& e){}
        try {
          double f = value.as<double>();
          dMap[p_name] = new double(f);
          continue;
        }
        catch (const std::exception& e){}
        try {
          std::string s = value.as<std::string>();
          sMap[p_name] = s;
          continue;
        }
        catch (const std::exception& e) {}
      } else { // Potentially this may go to minuit, if step is not zero.
        if (value.IsMap()) {
          yMap[p_name]=value;
        } else if (value.IsSequence() ) {
          size_t len = value.size();
	  vector<double> v(len);
	  try{
	    for (size_t i=0; i<len; i++)
	      v[i] = value[i].as<double>();
	    vMap[p_name] = v;
	  }catch(const YAML::TypedBadConversion<double>&ex){
	    try {
	      vector<string> v(len);
	      for (size_t i=0; i<len; i++)
		v[i] =  value[i].as<string>();
	      vsMap[p_name] = v;
	    }
	    catch(const YAML::TypedBadConversion<double>&ex){
	      cerr<<"[ERROR] parse_node_ failed to parse vector-parameter \""<<p_name<<"\":"<<value<<endl;
	      hf_errlog(18112100,"F: parse_node_ failed to parse sequence with non-double/string elements, see stderr");
	    };
	  }
        }
      }
    }
  }

  void createOutputDir(){
    string outputDir = "output";//default
    if(rootNode["OutputDirectory"]){
      outputDir = rootNode["OutputDirectory"].as<string>();
    }
    stringToFortran(coutdirname_.outdirname, 256, outputDir);

    //create the output directory
    //if that directory already exists, presumably from the previous run, rename it to output_OLD first
    //if output_OLD also exists, delete it
    if(fileExists(outputDir)){
      string oldOutputDir=outputDir+"_OLD";
      if(fileExists(oldOutputDir)){
        hf_errlog(1303201701, "W: Removing old results directory "+oldOutputDir);
        int ret = system(("rm -rf "+oldOutputDir).c_str());
      }
      hf_errlog(1303201702, "W: Backing up previous results to "+oldOutputDir);
      rename(outputDir.c_str(), oldOutputDir.c_str());
    }
    mkdir(outputDir.c_str(),0755);
  }

  void ParsToFortran(){

    // helper macros
#define FortAssignD(NameV,Struc)  if (gParameters.find(#NameV) != gParameters.end()) Struc.NameV = *gParameters[#NameV];
#define FortAssignS(NameV,Struc)  if (gParametersS.find(#NameV) != gParametersS.end()) strcpy(Struc.NameV,gParametersS[#NameV].c_str());
#define FortAssignI(NameV,Struc)  if (gParametersI.find(#NameV) != gParametersI.end()) Struc.NameV = gParametersI[#NameV];

    // CKM:
    FortAssignD(Vud,ckm_matrix_)
    FortAssignD(Vus,ckm_matrix_)
    FortAssignD(Vub,ckm_matrix_)
    FortAssignD(Vcd,ckm_matrix_)
    FortAssignD(Vcs,ckm_matrix_)
    FortAssignD(Vcb,ckm_matrix_)
    FortAssignD(Vtd,ckm_matrix_)
    FortAssignD(Vts,ckm_matrix_)
    FortAssignD(Vtb,ckm_matrix_)

    // Boson masses
    FortAssignD(Mz,boson_masses_)
    FortAssignD(Mw,boson_masses_)
    FortAssignD(Mh,boson_masses_)

    // Width
    FortAssignD(Wz,widths_)
    FortAssignD(Ww,widths_)
    FortAssignD(Wh,widths_)
    FortAssignD(Wtp,widths_)

    // EW couplings
    FortAssignD(Alphaem,ew_couplings_)

    //Weinberg angle
    if( rootNode["cos2thW"] or !rootNode["sin2thW"] ) {
      cerr<<"[ERROR] Weinberg angle must be given as sin2thW, not as cos2thW. Make sure sin2thW is specified in YAML steering and cos2thW is not."<<endl;
      hf_errlog(19090500, "F: Bad Weinberg angle: need sin2thW, not cos2thW, see stderr");
    }
    ew_couplings_.sin2thW = rootNode["sin2thW"].as<double>();
    ew_couplings_.cos2thW = 1. - ew_couplings_.sin2thW;

    // constants
    FortAssignD(Gf,constants_)
    FortAssignD(ConvFac,constants_)

    // OZ 26.04.18 some lines above do not do what expected because it is gf, convFac and alphaem in parameters.yaml, not Gf, ConvFac and Alphaem
    // as a result CC DIS cross section and integrated NC and CC cross sections are always zero with old interface
    // temporary fix: set these parameters manually
    // (maybe some other parameters are not assigned as well)
    if(gParameters.find("alphaem") != gParameters.end())
      ew_couplings_.Alphaem = *gParameters["alphaem"];
    if(gParameters.find("gf") != gParameters.end())
      constants_.Gf = *gParameters["gf"];
    if(gParameters.find("convFac") != gParameters.end())
      constants_.ConvFac = *gParameters["convFac"];

    //Fermion masses:
    FortAssignD(men,fermion_masses_)  // electron neutrino
    FortAssignD(mel,fermion_masses_)
    FortAssignD(mmn,fermion_masses_)
    FortAssignD(mmo,fermion_masses_)
    FortAssignD(mtn,fermion_masses_)
    FortAssignD(mta,fermion_masses_)
    FortAssignD(mup,fermion_masses_)
    FortAssignD(mdn,fermion_masses_)
    FortAssignD(mch,fermion_masses_)
    FortAssignD(mst,fermion_masses_)
    FortAssignD(mtp,fermion_masses_)
    FortAssignD(mbt,fermion_masses_)
  }

  const double*createConstantParameter(const string&n,double val){
    if(gParameters.count(n)>0){
      cerr<<"[ERROR] Redefinition of parameter \""<<n<<"\" with a constant parameter"<<endl;
      hf_errlog(19042101,"F: Parameter redefinition as a constant, see stderr");
    }
    double*p=new double;
    *p=val;
    //TODO: remember these created parameters and delete them
    //needs to be synchronized with how minimization parameters are deleted
    gParameters[n]=p;
    return p;
  }
  void createParameters(){
    using namespace std;
    try{
      YAML::Node parsNode=XFITTER_PARS::rootNode["Parameters"];
      if(!parsNode)return;//if there is no "Parameters:" node, just do not create any parameters
      if(!parsNode.IsMap()){
        hf_errlog(18091710,"F: Failed to create parameters: bad \"Parameters\" YAML node");
      }
      //Process parameters in alphabetical order
      //This is needed for MINUIT-assigned parameter IDs to be deterministic
      vector<string> parameterNames;
      for (const auto it:parsNode){
        try{
          parameterNames.push_back(it.first.as<string>());
        }catch(YAML::TypedBadConversion<string>){
          cerr<<"[ERROR] Failed to convert parameter name to string; trying to print it"<<endl;
          cerr<<it.first<<endl;
          hf_errlog(19063010, "F: Failed to convert parameter name to string, see stderr");
        }
      }
      sort(parameterNames.begin(),parameterNames.end());

      xfitter::BaseMinimizer*minimizer=xfitter::get_minimizer();
      vector<xfitter::DependentParameter>dependentParameters;

      for (auto it:parameterNames){
        string parameterName=it;
        stripString(parameterName);
        if(XFITTER_PARS::gParameters.find(parameterName)!=XFITTER_PARS::gParameters.end()){
          cerr<<"[ERROR] Redefinition of parameter \""<<parameterName<<"\""<<endl;
          hf_errlog(18112810,"F: Parameter redefinition, see stderr");
        }
        if(XFITTER_PARS::rootNode[parameterName].IsDefined()){
          cerr<<"[WARN] Redefinition of minimization parameter \""<<parameterName<<"\" at root of YAML steering. The value under \"Parameters\" will be used, and the one at global scope will be ignored"<<endl;
          hf_errlog(19042100,"W: Parameter redefinition at global scope, see stderr");
        }
        double value=nan("");
        double step=nan("");
        double min=nan("");
        double max=nan("");
        double pr_mean=nan("");
        double pr_sigma=nan("");
        YAML::Node pNode=parsNode[it];
        switch(pNode.Type()){
          case YAML::NodeType::Scalar:{
            string definition=pNode.as<string>();
            stripString(definition);
            if(definition=="DEPENDENT"||definition=="SUMRULE"){//This means that this parameter will be calculated using sum rules
              value=1;
              step=0;
            }else if(definition[0]=='='){//This means that dependence of parameter is given explicitly, e.g. "=2*B+1"
              value=1;
              step=0;
              dependentParameters.push_back({parameterName,definition});
            }else{//constant parameter defined by its value
              try{
                value=pNode.as<double>();
                step=0;
                break;
              }catch(YAML::TypedBadConversion<double>&ex){}
              cerr<<"[ERROR] Unable to parse definition of parameter \""<<parameterName<<"\": \""<<definition<<"\""<<endl;
              hf_errlog(18091712,"F: Bad parameter definition, see stderr");
            }
            break;}
          case YAML::NodeType::Sequence:{// interpret sequence as [value,step,min,max,priorVal,priorUnc]
            int size=pNode.size();
            if(size>0)value=pNode[0].as<double>();
            else break;
            if(size>1)step=pNode[1].as<double>();
            else break;
            if(size>2)min=pNode[2].as<double>();
            else break;
            if(size>3)max=pNode[3].as<double>();
            else break;
            if(size>4)pr_mean=pNode[4].as<double>();
            else break;
            if(size>5)pr_sigma=pNode[5].as<double>();
            else break;
            if(size>6){
              cerr<<"[ERROR] Too many arguments("<<size<<") in definition of parameter "<<parameterName<<"; expected at most 6"<<endl;
              hf_errlog(18091716,"F: Too many arguments in parameter definition, see stderr");
            }
            break;}
          case YAML::NodeType::Map:{
            for(YAML::const_iterator it=pNode.begin();it!=pNode.end();++it){
              string key=it->first.as<string>();
              if     (key=="value")   value   =it->second.as<double>();
              else if(key=="step")    step    =it->second.as<double>();
              else if(key=="min")     min     =it->second.as<double>();
              else if(key=="max")     max     =it->second.as<double>();
              else if(key=="pr_mean") pr_mean =it->second.as<double>();
              else if(key=="pr_sigma")pr_sigma=it->second.as<double>();
              else{
                cerr<<"[ERROR] Unknown key "<<key<<" in definition of parameter "<<parameterName<<endl;
                hf_errlog(18091711,"F: Unknown key in parameter definition, see stderr");
              }
            }
            break;
          }
          default:
            hf_errlog(18091710,"F: Failed to parse parameters definition node, expected Sequence or Map");
        }
        if(std::isnan(value)){
          value=1;
          hf_errlog(18091713,"W: Initial value not given for a parameter, assuming 1.0");
        }
        if(std::isnan(step)||step<0){
          step=fabs(value)*0.01;
          hf_errlog(18091714,"W: Step not given for a parameter, assuming 1%");
        }
        if(std::isnan(min)^std::isnan(max)){
          cerr<<"[ERROR] Only one of two bounds ("<<(std::isnan(min)?"upper":"lower")<<") is given for parameter "<<parameterName<<endl;
          hf_errlog(18091715,"F: Only one of two bounds is given for one of parameters, see stderr");
        }
        if(std::isnan(pr_mean)^std::isnan(pr_sigma)){
          if(std::isnan(pr_sigma)){
            pr_sigma=1;
            hf_errlog(18091716,"W: Prior mean value, but not sigma given for a parameter, assuming sigma=1.0");
          }else{
            cerr<<"[ERROR] Prior sigma value, but not prior mean fiven for parameter "<<parameterName<<endl;
            hf_errlog(18091717,"F: Prior sigma value, but not mean given for a parameter, see stderr");
          }
        }
        //Note that if we allocate memory here, it will never be freed
        double*bounds=nullptr;
        double*priors=nullptr;
        if(!std::isnan(min)){
          bounds=new double[2];
          bounds[0]=min;
          bounds[1]=max;
        }
        if(!std::isnan(pr_mean)){
          priors=new double[2];
          priors[0]=pr_mean;
          priors[1]=pr_sigma;
        }
        minimizer->addParameter(value,parameterName,step,bounds,priors);
      }
      xfitter::registerDependentParameters(dependentParameters);
    }catch(YAML::Exception&ex){
      cerr<<"[ERROR] YAML exception:\n"<<ex.what()<<"\n[/ERROR]"<<endl;
      hf_errlog(18091713,"F: YAML exception while creating parameters, details written to stderr");
    }
  }

  void createParameterisations(){
    //Read YAML steering and create the parameterisations that are defined there
    YAML::Node paramsNode=rootNode["Parameterisations"];
    if(!paramsNode)return;
    for(const auto it:paramsNode){
      string name=it.first.as<string>();
      try{
        YAML::Node definition=it.second;
        string classname=definition["class"].as<string>();
        gParameterisations[name] = (BasePdfParam*)createDynamicObject("pdfparam", classname, name);
      }catch(const YAML::InvalidNode&ex){
        cerr<<"[ERROR] Failed to create parameterisation \""<<name<<"\": bad YAML definition"<<endl;
        hf_errlog(20021500,"F: Bad parameterisation definition, see stderr");
      }
    }
    //Initialize
    for(const auto&pdfparam:gParameterisations){
      try{
        pdfparam.second->atStart();
      }catch(...){
        string name=pdfparam.first;
        cerr<<"[ERROR] Unhandled exception during initialization of parameterisation \""<<name<<"\". Check that its YAML definition is correct. Rethrowing the exception..."<<endl;
        throw;
      }
    }
  }
  void createDecompositions(){
    //Read YAML steering and create the decompositions that are defined there
    YAML::Node decompsNode=rootNode["Decompositions"];
    if(!decompsNode)return;
    for(const auto it:decompsNode){
      string name=it.first.as<string>();
      try{
        YAML::Node definition=it.second;
        string classname=definition["class"].as<string>();
        gPdfDecompositions[name] = (BasePdfDecomposition*)createDynamicObject("pdfdecomp", classname, name);
      }catch(const YAML::InvalidNode&ex){
        cerr<<"[ERROR] Failed to create decomposition \""<<name<<"\": bad YAML definition"<<endl;
        hf_errlog(20021501,"F: Bad decomposition definition, see stderr");
      }
    }
    //Initialize
    for(const auto&decomp:gPdfDecompositions){
      try{
        decomp.second->atStart();
      }catch(...){
        string name=decomp.first;
        cerr<<"[ERROR] Unhandled exception during initialization of decomposition \""<<name<<"\". Check that its YAML definition is correct. Rethrowing the exception..."<<endl;
        throw;
      }
    }
  }
  void createEvolutions(){
    //Read YAML steering and create the evolutions that are defined there
    YAML::Node evolsNode=rootNode["Evolutions"];
    if(!evolsNode)return;
    for(const auto it:evolsNode){
      string name=it.first.as<string>();
      try{
        YAML::Node definition=it.second;
        string classname = definition["class"].as<string>();
        gEvolutions[name] = (BaseEvolution*)createDynamicObject("evolution", classname, name);
      }catch(const YAML::InvalidNode&ex){
        cerr<<"[ERROR] Failed to create evolution \""<<name<<"\": bad YAML definition"<<endl;
        hf_errlog(20021502,"F: Bad evolution definition, see stderr");
      }
    }
    //Initialize
    for(const auto&evolution:gEvolutions){
      try{
        evolution.second->atStart();
      }catch(...){
        string name=evolution.first;
        cerr<<"[ERROR] Unhandled exception during initialization of evolution \""<<name<<"\". Check that its YAML definition is correct. Rethrowing the exception..."<<endl;
        throw;
      }
    }
  }

  std::string getParamFromNodeS(const std::string& name, const YAML::Node& node)
  {
    if(node[name].IsDefined())
      return node[name].as<string>();
    else {
      hf_errlog(19052019, "F: Undefined parameter \"" + name + "\" in node \"" + "\" requested as string");
      return "";
    }
  }
}

namespace xfitter{

BaseMinimizer* get_minimizer() {
  bool HasMinimizer  = ( XFITTER_PARS::gParametersS.find("Minimizer" ) != XFITTER_PARS::gParametersS.end() );
  bool HasMinimizers = ( XFITTER_PARS::gParametersVS.find("Minimizers" ) != XFITTER_PARS::gParametersVS.end() );

  if ( HasMinimizers && HasMinimizer ) {
    hf_errlog((int) 2203060601, "F: Both Minimizer and Minimizers present in parameters.yaml. Keep only one");
  }

  std::string name("");

  if (HasMinimizer)
    name = XFITTER_PARS::getParamS("Minimizer");
  else
    if ( XFITTER_PARS::gParametersS.find("__currentMinimizer" ) != XFITTER_PARS::gParametersS.end() )
      name = XFITTER_PARS::getParamS("__currentMinimizer" );
    else{
      name = XFITTER_PARS::gParametersVS.at("Minimizers")[0];
      XFITTER_PARS::gParametersS["__currentMinimizer"] = name;
    }

  // Check if already present
  if ( XFITTER_PARS::gMinimizer && XFITTER_PARS::gMinimizer->getName() == name ) {
    return  XFITTER_PARS::gMinimizer;  //already loaded
  }

  // copy original
  auto previousMinimizer =  XFITTER_PARS::gMinimizer;

  // else load, initialize and return
  XFITTER_PARS::gMinimizer =(BaseMinimizer*) createDynamicObject("minimizer", name);
  if (previousMinimizer)
    XFITTER_PARS::gMinimizer->CopyStateFromMinimizer(previousMinimizer);
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

void ensureMapValidity(const string&nodeName){
  //Report an error if a YAML map has duplicate keys
  //This is used for checking redefinition of parameterisations etc
  YAML::Node node=XFITTER_PARS::rootNode[nodeName];
  if(!node)return;
  if(!node.IsMap()){
    cerr<<"[ERROR] Node \""<<nodeName<<"\" is not a map"<<endl;
    hf_errlog(19040131,"F: Bad map node, see stderr");
  }
  set<string>keys;
  for(YAML::const_iterator it=begin(node);it!=end(node);++it){
    try{
      string key=it->first.as<string>();
      if(keys.count(key)!=0){
        cerr<<"[ERROR] Duplicate key \""<<key<<"\" in map \""<<nodeName<<'\"'<<endl;
        hf_errlog(19040132,"F: Duplicate key in a map node, see stderr");
      }
      keys.insert(key);
    }catch(YAML::TypedBadConversion<string>&ex){
      cerr<<"[ERROR] In map \""<<nodeName<<"\": failed to convert the following key to string:\n"<<it->first<<endl;
      hf_errlog(19040131,"F: Bad key in a map node, see stderr");
    }
  }
}

void parse_params_(){
  using namespace XFITTER_PARS;
  rootNode=loadYamlFile("parameters.yaml");
  expandIncludes(rootNode);
  createOutputDir();
  ensureMapValidity("Parameterisations");
  ensureMapValidity("Decompositions");
  ensureMapValidity("Evolutions");
  ensureMapValidity("byReaction");
  parse_node(rootNode,gParameters,gParametersI,gParametersS,gParametersV,gParametersVS,gParametersY);
  ParsToFortran();
  createParameters();
  createParameterisations();
  createDecompositions();
  createEvolutions();
}

// Store parameter to the map, fortran interface. Note that ref to the map travels from c++ to fortran and back:
void add_to_param_map_(map<std::string,double*> *map, double &value, int& global, char *name, int len) {
  string nam = name;
  const auto pos = nam.find(" ");
  if (pos < nam.size()) {
    nam.erase(pos);
  }

  if ( global>0 ) {
    XFITTER_PARS::gParameters[nam] = &value;
  }
  else {
    (*map)[nam] = &value;
  }
}

double getparamd_(const char* name,int len){
  char buff[128];
  memcpy( buff, &name[0], len);
  buff[len] = '\0';
  std::string key(buff);
  if (XFITTER_PARS::gParameters.find(key) != XFITTER_PARS::gParameters.end()) {
    return *XFITTER_PARS::gParameters[key];
  }
  else {
    return 0;
  }
}

int getparami_(const char* name,int len){
  char buff[128];
  memcpy( buff, &name[0], len);
  buff[len] = '\0';
  std::string key(buff);
  if (XFITTER_PARS::gParametersI.find(key) != XFITTER_PARS::gParametersI.end()) {
    return XFITTER_PARS::gParametersI[key];
  }
  else {
    return 0;
  }
}

void update_pars_fortran_() {
  XFITTER_PARS::ParsToFortran();
}
