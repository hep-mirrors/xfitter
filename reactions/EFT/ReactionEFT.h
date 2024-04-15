#ifndef xFitter_ReactionEFT
#define xFitter_ReactionEFT

#pragma once

#include "TermData.h"
#include "ReactionTheory.h"
#include "EFTTerm.h"
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

     @version 0.1
     @date 2023-05
  */

// class EFTReaction : public EFTTerm {
// public:
//  EFTReaction(vector<string> fname_list, ReactionTheory* reaction) : EFTTerm(fname_list, reaction) {};
// protected:
//   EFTReaction(vector<string> fname_list) : EFTTerm(fname_list) {}; // not public!
// };


class ReactionEFT : public ReactionTheory {
public:
  ReactionEFT(){};

public:
  std::map<int, EFTTerm* > EFT_terms;
  //
  virtual string getReactionName() const override { return  "EFT" ;};
  virtual void initTerm(TermData* td) override final;
  virtual void freeTerm(TermData*) override final;
  virtual void atStart() override {};
  virtual void compute(TermData*, valarray<double> &val, map<string, valarray<double> > &err) override;
protected:
  virtual int parseOptions(){ return 0;};
private:
  long int num_comp = 0; // number of iterations; used for debugging and profiling
  double time_comp = 0.0;  
};

#endif
