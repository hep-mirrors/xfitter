/*!
 @file xfitter_pars.cc
 @date Sun 16 April 2017
 @author SG

 Contains functions to read parameters.yaml, 
 global maps to store parameters,  and fortran interface functions.
 */

#include "xfitter_pars.h"
#include "xfitter_cpp.h"
#include <fstream>
#include <string.h>


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
}


namespace XFITTER_PARS {

  // Global vars:
  map <string, double*> gParameters;
  map <string, int>    gParametersI;
  map <string, string> gParametersS;
  map <string, vector<double> > gParametersV; ///< Vectors of double parameters
  map <string, YAML::Node > gParametersY;      ///< Store complete nodes for complex cases


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
	  if (value["step"]) {  
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
	    // no step, store as it is as a yaml node:
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
}

void parse_params_(){
  XFITTER_PARS::parse_file("parameters.yaml");
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
