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

  using namespace xfitter;
  xfitter::InitialPDFfunction getInputFunctionFromYaml(const YAML::Node&rootNode){
    YAML::Node node=rootNode["decomposition"];
    string name;
    try{
      name=node.as<string>();
    }catch(YAML::TypedBadConversion<string>ex){
      ostringstream s;s<<"W: YAML exception: "<<ex.what()<<"; while trying to extract decomposition name from node: "<<node<<"; using default decomposition name";
      hf_errlog(18082930,s.str()); 
      name=getDefaultDecompositionName();
    }
    return xfitter::get_pdfDecomposition(name)->f0();
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
  bool is_file_exist(const char *fileName)
  {
    std::ifstream infile(fileName);
    return infile.good();
  }

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
  


  // Parse @param node and return maps
  void parse_node(const YAML::Node& node, 
                  std::map<string,double*>& dMap, 
                  std::map<string,int>& iMap, 
                  std::map<string,string>& sMap, 
                  std::map<string,vector<double> >& vMap,
                  std::map<string,YAML::Node> & yMap )
  {
    for ( YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
      YAML::Node key = it->first;
      YAML::Node value = it->second;
      
      // parameter name
      string p_name = key.as<string>();
    
      //  Check if asked to read another file:
      if ( p_name == "include" ) {
  auto fileName =  value.as<string>();
  if (is_file_exist(fileName.c_str())) {
    parse_file( fileName );
  }
  else {
    // Now try default location:
    if (is_file_exist((PREFIX +string("/")+fileName).c_str())) {
      parse_file( PREFIX +string("/")+fileName );
    }
    else {
      string msg = "F: Include Yaml parameters file "+fileName+" not found";
      hf_errlog_(17041601,msg.c_str(), msg.size());
    }
  }
      }
      
      if (value.IsScalar()) {
  // Alright, store directly
  // Try to read as int, float, string:
  try {
    int i = value.as<int>();
    iMap[p_name] = i;
    continue;
  }
  catch (const std::exception& e) {
  }
  
  try {
    double f = value.as<double>();
    dMap[p_name] = new double(f);
    continue;
  }
  catch (const std::exception& e) {
  }
  
  try {
    std::string s = value.as<std::string>();
    sMap[p_name] = s;
    continue;
  }
  catch (const std::exception& e) {
  }      
      }
      else { // Potentially this may go to minuit, if step is not zero.
  
  if (value.IsMap()) {

    // Check if this is a minimisation block, true if step is present
    if (value["step"] || value["value"]) {  
      // Defaults
      double val = 0;
      double step = 0;
      double minv  = 0;
      double maxv  = 0;
      double priorVal = 0;
      double priorUnc = 0;
      int add = true;
      
      if (value["value"]) {
        val = value["value"].as<double>();
      }
      else {
        string text = "F: missing value field for parameter " + p_name;
        hf_errlog_(17032401,text.c_str(),text.size());       
      }
      if (value["step"]) step = value["step"].as<double>();
      if (value["prior"]) priorVal = value["prior"].as<double>();
      if (value["priorUnc"]) priorUnc = value["priorUnc"].as<double>();
      if (value["min"]) minv = value["min"].as<double>();
      if (value["max"]) maxv = value["max"].as<double>();
      // Goes to fortran
      addexternalparam_(p_name.c_str(),  val, step, minv, maxv,
            priorVal, priorUnc, add, &dMap, p_name.size());
    }
    else {
      // no step or value, store as it is as a yaml node:
      yMap[p_name] = value;
    }
  }
  
  else if (value.IsSequence() ) {
    size_t len = value.size();
    vector<double> v(len);
    for (size_t i=0; i<len; i++) {
      v[i] = value[i].as<double>();
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

  const std::function<void(double const& x, double const& Q, double* pdfs)>  retrieveXfxQArray(const std::string& name) {
    return gXfxQArrays.at(name);
  }
  void createParameters(){
    using namespace std;
    try{
      YAML::Node parsNode=rootNode["Parameters"];
      if(!parsNode.IsMap()){
        hf_errlog(18091710,"F: Failed to create parameters: bad \"Parameters\" YAML node");
      }
      xfitter::BaseMinimizer*minimizer=xfitter::get_minimizer();
      for(YAML::const_iterator it=parsNode.begin();it!=parsNode.end();++it){
        string parameterName=it->first.as<string>();
        double value=nan("");
        double step=nan("");
        double min=nan("");
        double max=nan("");
        double pr_mean=nan("");
        double pr_sigma=nan("");
        //TODO: bounds and priors
        YAML::Node pNode=it->second;
        switch(pNode.Type()){
          case YAML::NodeType::Scalar:{//Should be a special string DEPENDENT
            string key=pNode.as<string>();
            if(key=="DEPENDENT"){
              value=1;
              step=0;//This means that this parameter will be calculated using sum rules
            }
            else hf_errlog(18091712,"F: Bad parameter definition");
            break;}
          case YAML::NodeType::Sequence:{// interpret sequence as [value,step,min,max]
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
              if   (key=="value")value=it->second.as<double>();
              else if(key=="step")step=it->second.as<double>();
              else if(key=="min")min=it->second.as<double>();
              else if(key=="max")max=it->second.as<double>();
              else if(key=="pr_mean")pr_mean=it->second.as<double>();
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
    }catch(YAML::Exception&ex){
      cerr<<"[ERROR] YAML exception:\n"<<ex.what()<<"\n[/ERROR]"<<endl;
      hf_errlog(18091713,"F: YAML exception while creating parameters, details written to stderr");
    }
  }
}

void parse_params_(){
  XFITTER_PARS::parse_file("parameters.yaml");
  XFITTER_PARS::createParameters();
  XFITTER_PARS::ParsToFortran();
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
