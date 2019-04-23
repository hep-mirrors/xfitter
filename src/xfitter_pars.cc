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
#include <fstream>
#include <string.h>
#include <cmath>
#include "BaseEvolution.h"
#include "BasePdfDecomposition.h"
#include "BaseMinimizer.h"
#include"dependent_pars.h"

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
  // Update of EWK/QCD parameters, can be fitted at each iteration.
  void update_pars_fortran_(); 
}


namespace XFITTER_PARS {

  // Global vars:
  xfitter::BaseMinimizer* gMinimizer(nullptr);
  YAML::Node rootNode;
  map <string, double*> gParameters;
  map <string, int>    gParametersI;
  map <string, string> gParametersS;
  map <string, vector<double> > gParametersV; ///< Vectors of double parameters
  map <string, YAML::Node > gParametersY;      ///< Store complete nodes for complex cases
  
  map<string,std::function<void(double const& x, double const& Q, double* pdfs)> > gXfxQArrays;

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
      cerr<<"[ERROR] Double parameter \""<<name<<"\" does not exist; rethrowing out_of_range"<<endl;
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
    string name;
    try{
      name=node.as<string>();
    }catch(YAML::TypedBadConversion<string>ex){
      ostringstream s;s<<"W: YAML exception: "<<ex.what()<<"; while trying to extract decomposition name from node: "<<node<<"; using default decomposition name";
      hf_errlog(18082930,s.str()); 
      name=getDefaultDecompositionName();
    }
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

  std::string getParameterS(std::string name) {
    auto search = gParametersS.find(name);
    if ( search != gParametersS.end() ) {
      return search->second;
    }
    else {
      hf_errlog(18071301,"W: string parameter "+name+" not found"); 
      return ""; // not found
    }
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
// Helper function
bool fileExists(const string&fileName){
  return std::ifstream(fileName).good();
}

/*
  void parse_file(const std::string& name)
  {
    try {
      if ( ! std::ifstream(name).good()) {
        string text = "F: Problems opening parameters file " + name;
        hf_errlog_(18032001,text.c_str(), text.size());
      }
      YAML::Node node = YAML::LoadFile(name);
      parse_node(node, gParameters, gParametersI, gParametersS, gParametersV, gParametersY);
      //HACKY way to get rootNode, pending include rewrite
      if(rootNode.IsNull())rootNode=node;
    }
    catch (const std::exception& e) {
      std::cout << e.what() << std::endl;
      string text = "F: Problems reading yaml-parameters file: " + name + " : "+e.what();
      hf_errlog_(17032503,text.c_str(),text.size());
    }
    
    return;
  }
*/
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
void expandIncludes(YAML::Node&node,unsigned int recursionLimit=256){
  if(recursionLimit==0){
    hf_errlog(18092605,"F: Recursion limit reached while handling includes");
  }
  if(!node.IsMap())return;//maybe this should even be an error
  vector<YAML::Node>include_keys;
  for(YAML::iterator it=node.begin();it!=node.end();++it){
    YAML::Node&key=it->first;
    YAML::Node&val=it->second;
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
        string prefix_filename=PREFIX+string("/yaml/")+filename;
        if(fileExists(prefix_filename))filename=prefix_filename;
        else file_not_found=true;
      }
      if(file_not_found){
        cerr<<"[ERROR] YAML include file "<<filename<<" not found"<<endl;
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
          clog<<"[INFO] Option "<<key<<"="<<it->second<<" included from file "<<filename<<" is overridden by locally defined option "<<key<<"="<<node[key]<<endl;
          hf_errlog(18092604,"I: locally defined setting overrides included, see stdlog");
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
            for (size_t i=0; i<len; i++) {
              v[i] = value[i].as<double>();
            }
          }catch(const YAML::TypedBadConversion<double>&ex){
            cerr<<"[ERROR] parse_node_ failed to parse vector-parameter \""<<p_name<<"\":"<<value<<endl;
            hf_errlog(18112100,"F: parse_node_ failed to parse sequence with non-double elements, see stderr");
          }
          vMap[p_name] = v;
        }
      }
    }
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
    FortAssignD(sin2thW,ew_couplings_)
    FortAssignD(cos2thW,ew_couplings_)

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

    // Steering
    FortAssignS(hf_scheme,steering_)

  }

  void registerXfxQArray(const string& name, std::function<void(double const& x, double const& Q, double* pdfs)>  xfxArray) {
    gXfxQArrays[name] = xfxArray;
  }
  //TODO: delete this:
  const std::function<void(double const& x, double const& Q, double* pdfs)>  retrieveXfxQArray(const std::string& name) {
    return gXfxQArrays.at(name);
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
      if(!parsNode.IsMap()){
        hf_errlog(18091710,"F: Failed to create parameters: bad \"Parameters\" YAML node");
      }
      xfitter::BaseMinimizer*minimizer=xfitter::get_minimizer();
      vector<xfitter::DependentParameter>dependentParameters;
      for(YAML::const_iterator it=parsNode.begin();it!=parsNode.end();++it){
        string parameterName=it->first.as<string>();
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
        YAML::Node pNode=it->second;
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
            break;}
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
}
void ensureMapValidity(const string&nodeName){
  //Report an error if a YAML map has duplicate keys
  //This is used for checking redefinition of parameterisations etc
  YAML::Node node=XFITTER_PARS::rootNode[nodeName];
  if(!node){
    cerr<<"[ERROR] Necessary node \""<<nodeName<<"\" does not exist"<<endl;
    hf_errlog(19040134,"F: Necessary node does not exist, see stderr");
  }
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
  ensureMapValidity("Parameterisations");
  ensureMapValidity("Decompositions");
  ensureMapValidity("Evolutions");
  parse_node(rootNode,gParameters,gParametersI,gParametersS,gParametersV,gParametersY);
  createParameters();
  ParsToFortran();
}

// Store parameter to the map, fortran interface. Note that ref to the map travels from c++ to fortran and back:
void add_to_param_map_(map<std::string,double*> *map, double &value, int& global, char *name, int len) {
  string nam = name;
  nam.erase(nam.find(" "));

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

void update_pars_fortran_() {
  XFITTER_PARS::ParsToFortran();
}
