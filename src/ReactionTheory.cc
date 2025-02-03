/*!
 @file ReactionTheory.cc
 @date Thu Jan 21 2016
 @author Andrey Sapronov <sapronov@ifh.de>

 Contains implementations of ReactionTheory class member functions.
 */
#include "ReactionTheory.h"
#include <unistd.h>
#include "xfitter_cpp.h"
void ReactionTheory::atStart(){}
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
  _ncpu = td->_ncpu;
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
