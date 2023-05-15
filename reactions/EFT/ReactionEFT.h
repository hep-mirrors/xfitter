#ifndef xFitter_ReactionEFT
#define xFitter_ReactionEFT

#pragma once

#include "ReactionTheory.h"
#include "EFTReader.h"
//#include <memory> 
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>

  /**
     @class' ReactionEFT

     @brief A wrapper class for EFT reaction 

     Based on the ReactionTheory class
     Adapted from ReactionCIJET

     @version 0.3
     @date 2016-12-06
  */

  // a7: why we need EFTReaction when we already have EFTReader?
class EFTReaction : public EFTReader {
public:
  EFTReaction(string name, ReactionTheory* reaction) : EFTReader(name, reaction) {};
protected:
  EFTReaction(string name) : EFTReader(name) {}; // not public!
};


class ReactionEFT : public ReactionTheory {
public:
  ReactionEFT(){};

public:
  virtual string getReactionName() const override { return  "EFT" ;};
  virtual void initTerm(TermData* td) override final;
  virtual void atStart() override {};
  // virtual void freeTerm(TermData*) override final; // delete the pineappl grids.
  virtual void compute(TermData*, valarray<double> &val, map<string, valarray<double> > &err) override;
protected:
  virtual int parseOptions(){ return 0;};
  
  // std::map<int,std::vector<EFTReaction*> > EFT_all_dataset;
  std::map<int, EFTReaction* > EFT_all_dataset;
  vector<string> name_EFT_param; 
};

#endif
