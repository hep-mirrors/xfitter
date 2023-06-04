#include "xfitter_steer.h"
#include "xfitter_pars.h"
#include "xfitter_cpp_base.h"
#include"xfitter_cpp.h"
#include <iostream>
#include<fstream>
#include <yaml-cpp/yaml.h>
#include <Profiler.h>
#include"BasePdfDecomposition.h"
#include"BaseEvolution.h"
#include"BaseMinimizer.h"
using std::string;
using std::cerr;

extern std::map<string,string> gReactionLibs;


namespace xfitter {
BasePdfParam*getParameterisation(const string&name){
  try{
    return XFITTER_PARS::gParameterisations.at(name);
  }catch(const std::out_of_range&){
    cerr<<"[ERROR] Parameterisation \""<<name<<"\" does not exist; rethrowing out_of_range"<<endl;
    throw;
  }
}
BasePdfDecomposition*get_pdfDecomposition(string name){
  if(name=="")name=XFITTER_PARS::getDefaultDecompositionName();
  try{
    return XFITTER_PARS::gPdfDecompositions.at(name);
  }catch(const std::out_of_range&){
    cerr<<"[ERROR] Decomposition \""<<name<<"\" does not exist; rethrowing out_of_range"<<endl;
    throw;
  }
}
BaseEvolution*get_evolution(string name){
  if(name=="")name=XFITTER_PARS::getDefaultEvolutionName();
  try{
    return XFITTER_PARS::gEvolutions.at(name);
  }catch(const std::out_of_range&){
    cerr<<"[ERROR] Evolution \""<<name<<"\" does not exist; rethrowing out_of_range"<<endl;
    throw;
  }
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
      hf_errlog(16042803,"E: Error matrix not accurate");
      break;
    case ConvergenceStatus::FORCED_POSITIVE:
      hf_errlog(16042804,"E: Error matrix forced positive");
      break;
    case ConvergenceStatus::SUCCESS:
      hf_errlog(16042802,"I: Fit converged");
      break;
    case ConvergenceStatus::NO_CONVERGENCE:
      hf_errlog(16042805,"E: No convergence");
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
