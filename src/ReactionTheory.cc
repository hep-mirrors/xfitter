/*!
 @file ReactionTheory.cc
 @date Thu Jan 21 2016
 @author Andrey Sapronov <sapronov@ifh.de>

 Contains implementations of ReactionTheory class member functions.
 */
#include "ReactionTheory.h"
#include <unistd.h>
#include "xfitter_cpp.h"
#include "xfitter_pars.h"
void ReactionTheory::atStart(){
  // todo: need to generalise this code, I want to get "threads" from reaction or global
  _ncpu = 1;
  const std::string BY_REACTION = "byReaction";
  const std::string parName = "threads";
  if(XFITTER_PARS::gParameters.count(parName) > 0 or XFITTER_PARS::rootNode[parName].IsDefined()) {
    _ncpu = XFITTER_PARS::rootNode[parName].as<int>();
  }
  YAML::Node byReactionNode = XFITTER_PARS::rootNode[BY_REACTION];
  if (byReactionNode.IsMap()) {
    YAML::Node reactionNode = byReactionNode[getReactionName()];
    if (reactionNode.IsMap() && reactionNode[parName].IsDefined()) {
      _ncpu = reactionNode[parName].as<int>();
    }
  }
  if (_ncpu == -1) {
      _ncpu = sysconf(_SC_NPROCESSORS_ONLN);
      hf_errlog(2023061401,"I: Will use "+std::to_string(_ncpu)+" threads");
  }
}
void ReactionTheory::atIteration(){}
void ReactionTheory::initTerm(TermData*td){
  // parallel
  td->_ncpu = 1;
  if (td->hasParam("threads"))
    td->_ncpu = td->getParamI("threads");
  if (td->_ncpu == -1) {
      td->_ncpu = sysconf(_SC_NPROCESSORS_ONLN);
      hf_errlog(2023061401,"I: Will use "+std::to_string(td->_ncpu)+" threads");
  }
  _task_distr = ForkPool::TaskDistribution::cyclic;
  if (td->hasParam("parallel_task_distribution")) {
    _task_distr = ForkPool::getTaskDistrFromStr(td->getParamS("parallel_task_distribution"));
  }
  if (td->hasParam("computeAtIteration")) {
    _flagComputeAtIteration = td->getParamI("computeAtIteration");
  }
};
void ReactionTheory::freeTerm(TermData*td){};
void ReactionTheory::reinitTerm(TermData*td){
  freeTerm(td);
  initTerm(td);
};
void ReactionTheory::atFCN3(){}
void ReactionTheory::atMakeErrorBands(int i){}
int ReactionTheory::getNCPU(TermData* td) {
  int ncpu = 1;
  if (td->hasParam("threads")) {
    ncpu = td->getParamI("threads");
  }
  if (ncpu == -1) {
    ncpu = sysconf(_SC_NPROCESSORS_ONLN);
    hf_errlog(2026061200,"I: Will use "+std::to_string(ncpu)+" threads");
  }
  return ncpu;
}
