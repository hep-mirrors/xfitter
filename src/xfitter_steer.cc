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
#include "BasePdfParam.h"
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif


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

  pid_t xf_fork(int NCPU) {
    auto out = fork();
    if (out == 0) {
      // Child: tear down inherited OpenMP pool. fork() kills all non-calling
      // threads, so libgomp's worker pool is zombie in the child - the next
      // parallel region would deadlock.
      // setenv is kept as an additional hint (some
      // runtimes re-read env vars). This is just a lazy safety net -
      // the real fix is 
      // TODO: to not fork() when OpenMP threads > 1 (see xf_ncpu).
#ifdef _OPENMP
      omp_set_num_threads(1);
#endif
      setenv("OMP_NUM_THREADS", "1", 1);
      setenv("OPENBLAS_NUM_THREADS", "1", 1);


      // Get currently availiable number of CPUs
      auto it=XFITTER_PARS::gParametersI.find("NCPUmax");
      if(it != XFITTER_PARS::gParametersI.end()){
	int nCPUmax = XFITTER_PARS::gParametersI.at("NCPUmax");
	if (nCPUmax>0) {
	  nCPUmax = nCPUmax / NCPU;
	  if (nCPUmax == 0) {
	    nCPUmax = 1;
	  }
	  XFITTER_PARS::gParametersI.at("NCPUmax") = nCPUmax;
	}
      }
    }
    return out;
  }
  
  const int xf_ncpu(int NCPU) {
    if (NCPU == 0) return 0;

    // If OpenMP threading is active, shared-memory parallelism in GetChisquare
    // already uses all requested cores. Forking children that each also try
    // to thread would oversubscribe. Force fork counts to 1 unless the user
    // opts out via OpenMP.allowForkWithThreads: true. 
    // TODO: allowForkWithThreads is untested for now.
    auto itOmp = XFITTER_PARS::gParametersI.find("OMP_NUM_THREADS");
    if (itOmp != XFITTER_PARS::gParametersI.end() && itOmp->second > 1 && NCPU > 1) {
      auto itAllow = XFITTER_PARS::gParametersI.find("allowForkWithThreads");
      bool allow = (itAllow != XFITTER_PARS::gParametersI.end() && itAllow->second != 0);
      if (!allow) {
        static bool warned = false;
        if (!warned) {
          hf_errlog(2024040103,
                    "I: OpenMP threads="+std::to_string(itOmp->second)+
                    " active; forcing fork-parallel counts to 1 "
                    "(set OpenMP.allowForkWithThreads: true to override).");
          warned = true;
        }
        return 1;
      }
    }

    auto it=XFITTER_PARS::gParametersI.find("NCPUmax");
    if(it != XFITTER_PARS::gParametersI.end()){

      int nCPUmax = XFITTER_PARS::gParametersI.at("NCPUmax");

      // determine automatically
      if (nCPUmax < 0) {
	nCPUmax = sysconf(_SC_NPROCESSORS_ONLN);
	hf_errlog(2023111601,"I: Will use "+std::to_string(nCPUmax)+" maximum nCPU");
	XFITTER_PARS::gParametersI.at("NCPUmax") = nCPUmax;
      }

      if (nCPUmax > 0) {
	return min(NCPU,nCPUmax);
      }
      else {
	return NCPU;
      }
    }
    else {
      return NCPU;
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

bool updateMinimizer() {
  if ( XFITTER_PARS::gParametersVS.find("Minimizers" ) == XFITTER_PARS::gParametersVS.end() )
    return false;
  auto currentMinimizer =std::find( XFITTER_PARS::gParametersVS["Minimizers"].begin(),
				    XFITTER_PARS::gParametersVS["Minimizers"].end(),
				    XFITTER_PARS::gParametersS["__currentMinimizer"]);

  if (std::next(currentMinimizer) == XFITTER_PARS::gParametersVS["Minimizers"].end())
    return false;
  else {
    XFITTER_PARS::gParametersS["__currentMinimizer"] = *std::next(currentMinimizer);
    // update it and init
    auto mini = xfitter::get_minimizer();
    // We need to re-initialize parameterisations (since parameters may moved)
    for(const auto& pdfparam : XFITTER_PARS::gParameterisations){
      pdfparam.second->atStart();
    }
    return true;
  }
}


void run_minimizer_() {
  /// get profiler too
  auto *prof = new xfitter::Profiler();

  prof->doProfiling();

  do {
    auto mini = xfitter::get_minimizer();
    mini->doMinimization();
  }
  while ( updateMinimizer() );
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

