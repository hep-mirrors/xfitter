#pragma once

#include <string>
#include <sys/types.h>

/**
     @class xfitter_steer

     @brief basic building blocks of xFitter running sequence

     @version 0.1
     @date 2018-07-13
*/

class ReactionTheory;

namespace xfitter
{
  class BaseEvolution;
  class BasePdfDecomposition;
  class BasePdfParam;
  class BaseMinimizer;

  ///Get evolution by its name
  //If name=="", returns default evolution
  BaseEvolution* get_evolution(std::string name="");  
  ///Get decomposition by its name
  BasePdfDecomposition* get_pdfDecomposition(std::string name="");
  BasePdfParam*getParameterisation(const std::string&name="");
  ///Get minimizer
  BaseMinimizer*get_minimizer();
  ///Get reaction
  ReactionTheory* getReaction(const std::string& name);


  //Call atConfiguration change for all evolutions and decompositions
  //Call this after changing some of configuration parameters
  //Used by Profiler
  void updateAtConfigurationChange();

  //When fortran code accesses pdfs, it accesses this default evolution
  //extern BaseEvolution*defaultEvolution;
  BaseEvolution* defaultEvolutionInstance();

  // Fork keeping track of present CPUs
  pid_t xf_fork(int NCPU); //perform a fork, reduce the pull of CPUs (changes gParametersI["NCPUmax"])
  const int xf_ncpu(int NCPU); // get number of CPUs which can be used, does not modify enviroment. 
}
