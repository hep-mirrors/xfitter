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
  const std::string BY_REACTION = "byReaction";
  const std::string parName = "threads";
  YAML::Node byReactionNode = XFITTER_PARS::rootNode[BY_REACTION];
  if (byReactionNode.IsMap()) {
    YAML::Node reactionNode = byReactionNode[getReactionName()];
    if (reactionNode.IsMap()) {
      _ncpu = reactionNode[parName].as<int>();
    }
  }
  else if(XFITTER_PARS::gParameters.count(parName) > 0 or XFITTER_PARS::rootNode[parName].IsDefined()) {
    _ncpu = XFITTER_PARS::rootNode[parName].as<int>();
  }
  else {
    _ncpu = 1;
  }
}
void ReactionTheory::atIteration(){}
void ReactionTheory::initTerm(TermData*td){
  // parallel
  td->_ncpu = 1;
  if (td->hasParam("threads"))
    td->_ncpu = td->getParamI("threads");
  printf("td->_ncpu = %d\n", td->_ncpu);
  if (td->_ncpu == -1) {
      td->_ncpu = sysconf(_SC_NPROCESSORS_ONLN);
      hf_errlog(2023061401,"I: Will use "+std::to_string(td->_ncpu)+" threads");
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
