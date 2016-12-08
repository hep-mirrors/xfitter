
#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactiontestZMVFNS

  @brief A wrapper class for testZMVFNS reaction 

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2016-12-08
  */

class ReactiontestZMVFNS : public ReactionTheory
{
  public:
    ReactiontestZMVFNS(){};

//    ~ReactiontestZMVFNS(){};
//    ~ReactiontestZMVFNS(const ReactiontestZMVFNS &){};
//    ReactiontestZMVFNS & operator =(const ReactionAtestZMVFNS &r){return *(new ReactiontestZMVFNS(r));};

  public:
    virtual string getReactionName() const { return  "testZMVFNS" ;};
    int initAtStart(const string &); 
    virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err);
  protected:
    virtual int parseOptions(){ return 0;};
};

