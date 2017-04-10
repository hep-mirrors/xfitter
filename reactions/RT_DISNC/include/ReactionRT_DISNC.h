
#pragma once

#include "ReactionBaseDISNC.h"

/**
   @class' ReactionRT_DISNC

   @brief A wrapper class for RT_DISNC reaction 

   Based on the ReactionTheory class. Reads options produces 3d cross section.
  
   @version 0.1
   @date 2017-04-10
*/

class ReactionRT_DISNC : public ReactionBaseDISNC
{
 private:
  typedef ReactionBaseDISNC Super;
 public:
  ReactionRT_DISNC(){};
 public:
  virtual string getReactionName() const { return  "RT_DISNC" ;};
  int initAtStart(const string &); 
  virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err);
  virtual void initAtIteration() override; 

 protected:
  virtual void F2 BASE_PARS override;
  virtual void F2gamma_RT BASE_PARS;
};

