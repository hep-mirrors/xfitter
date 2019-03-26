#pragma once
#include"ReactionTheory.h"
#include"Variant.h"
/*
TermData is a class that provides interface to all parameters for a given reaction term
TermData is responsible for:
*Providing interface to reaction parameters, resolving overwrites and global vs specific parameters
*Reading special parameters "evolution", "evolution1" and "evolution2", get pointer[s] to correct evolution[s]
*Switch XFXlike and other wrappers

One instance of TermData exists for each reaction term.
The class TheorEval manages all instances of TermData and passes them to ReactionTheory

It is intended to be used by ReactionTheory
*/
namespace XFITTER_PARS{
class TermData{
public:
  TermData(unsigned id,ReactionTheory*);//add any other parameters TheorEval needs to pass
  //Unique id of this term
  //In the past this was calculated as 1000*datasetID+number_of_term (see TheorEval::initReactionTerm)
  //In ReactionTheory called it was also incorrectly called dataSetID
  const unsigned id;
  //Get reaction parameter by its name, taking into account that term-specific parameters overshadow global etc.
  Variant getParam(string)const;
  ReactionTheory*reaction;
  void actualizeWrappers();//see wrappers below
  BaseEvolution*getEvolution()const;//returns first evolution for this term
  BaseEvolution*getEvolution(int i)const;//i is either 0 for evolution1, or 1 for evolution2
  //The following pointer can be used by ReactionTheory to store some additional data
  //for each reaction term. It should be managed by ReactioTheory only, do not touch it from elsewhere
  void*reactionData=nullptr;
  //Dataset*dataset()const; //for future, after we decide on a Dataset class
  //We should consider adding any other getters that could be useful to ReactionTheory
private:
  //TODO: Implement TermData
}
/* Wrappers
  If a ReacationTheory uses these wrapper, it must call TermData::actualizeWrappers
  for each term at each iteration, to make sure the wrappers wrap the correct evolution
*/
}
extern "C"{
void PDF_xfx_wrapper (const double&x,const double&Q,double*results);
void PDF_xfx_wrapper1(const double&x,const double&Q,double*results);
void  AlphaS_wrapper (const double&Q);
}
