#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsubobject-linkage"


/*
   @file Q0decompostion.cc
   @date 2025-04-25
   @author SG
*/

#include "Q0decomposition.h"
#include "xfitter_pars.h"
#include "xfitter_cpp_base.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

bool is_file_exist(const string& fileName)
{
  std::ifstream infile(fileName);
  return infile.good();
}

namespace xfitter {
  //For dynamic loading:
  extern "C" Q0decomposition*create(const char*name){
    return new Q0decomposition(name);
  }
  
  Q0decomposition::Q0decomposition(const char*name):BasePdfDecomposition{name}{}

  Q0decomposition::~Q0decomposition(){}

  const char*Q0decomposition::getClassName()const{return"Q0decomposition";}

  void Q0decomposition::atStart(){
    YAML::Node pars=XFITTER_PARS::getDecompositionNode(_name);
    string fileName;
    try{
      fileName=pars["fileName"].as<std::string>();
    }catch(YAML::TypedBadConversion<std::string>&ex){
      std::cerr<<"[ERROR] Bad file name given for Q0 decomposition "<<_name<<": set=\n"<<pars<<"\n[/ERROR]"<<std::endl;
      hf_errlog(24042501,"F: In Q0decomposition: failed to convert YAML node \"set\" to string, see stderr");
    }

    // check if exists first
    if (is_file_exist(fileName)) {
      const auto tab=YAML::LoadFile(fileName);

      // Check consistency of the evolution parameters:
      using XFITTER_PARS::getParamD;
      using XFITTER_PARS::getParamS;
      if (tab["Order"])
	if (tab["Order"].as<std::string>() != getParamS("Order"))
	  hf_errlog(24052601,"F:Q0decomposition: Inconsistent Order of evolution "
		    +tab["Order"].as<std::string>()+ " vs "+ getParamS("Order"));

      for (const std::string& par : std::array<std::string,6>{"Q0","Mz","alphas","mch","mbt","mtp"} ) {
	if (tab[par])
	  if (abs(tab[par].as<double>() - *getParamD(par))>1.0e-8)
	    hf_errlog(24052602,"F:Q0decomposition: Inconsistent evolution parameter value for "+par);
      }
      
      std::vector<double> xGrid;
      if (tab["xGrid"] && tab["xGrid"].IsSequence()) {
	std::cout << "ha \n";
	xGrid = tab["xGrid"].as<std::vector<double>>();
      }
      else {
	hf_errlog(24042503,"F: Q0decomposition: failed to find xGrid node"); 
      }
      // loop over all PDFs
      for (int i : std::array<int,13> {-6,-5,-4,-3, -2, -1, 1, 2, 3, 4, 5, 6, 21}) {
	std::vector<double> p;
	std::cout << "here " << i << std::endl;
	char mess[20];
	sprintf(mess,"%d",i);
	if (tab[mess] && tab[mess].IsSequence()) {
	  p = tab[mess].as<std::vector<double>>();
	  _PDFsAtQ0[i].set_points(xGrid,p);
	}
	else {
	  hf_errlog(24042504,"F: Q0decomposition: failed to find PDF "+ std::string(mess)); 
	}
	
      }
    }
    else {
      hf_errlog(24042501,"F: In Q0decomposition: failed to open file "+fileName);
    }
  }

  std::map<int,double>Q0decomposition::xfxMap(double x)const{
    std::map<int,double> out;
    for (const auto& entry : _PDFsAtQ0) {
      out[entry.first] = entry.second(x);
    }
    return out;
  }
}
#pragma GCC diagnostic pop
