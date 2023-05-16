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

// class EFTReaction : public EFTReader {
// public:
//  EFTReaction(vector<string> fname_list, ReactionTheory* reaction) : EFTReader(fname_list, reaction) {};
// protected:
//   EFTReaction(vector<string> fname_list) : EFTReader(fname_list) {}; // not public!
// };


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
  
  std::map<int, EFTReader* > EFT_terms;
  vector<string> name_EFT_param; 
  int debug = -1;
};

#endif
