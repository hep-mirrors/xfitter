/*!
 @file ReactionTheory.cc
 @date Thu Jan 21 2016
 @author Andrey Sapronov <sapronov@ifh.de>

 Contains implementations of ReactionTheory class member functions.
 */
#include "ReactionTheory.h"
void ReactionTheory::atStart(){}
void ReactionTheory::atIteration(){}
void ReactionTheory::initTerm(TermData*td){};
void ReactionTheory::freeTerm(TermData*td){};
void ReactionTheory::reinitTerm(TermData*td){
  freeTerm(td);
  initTerm(td);
};
void ReactionTheory::atFCN3(){}
void ReactionTheory::atMakeErrorBands(int i){}
