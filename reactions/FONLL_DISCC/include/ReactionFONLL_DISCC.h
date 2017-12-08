
#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionFONLL_DISCC

  @brief A wrapper class for FONLL_DISCC reaction 

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-11-29
  */

class ReactionFONLL_DISCC : public ReactionTheory
{
  public:
    ReactionFONLL_DISCC(){};

//    ~ReactionFONLL_DISCC(){};
//    ~ReactionFONLL_DISCC(const ReactionFONLL_DISCC &){};
//    ReactionFONLL_DISCC & operator =(const ReactionAFONLL_DISCC &r){return *(new ReactionFONLL_DISCC(r));};

  public:
    virtual string getReactionName() const { return  "FONLL_DISCC" ;};
    int initAtStart(const string &); 
    virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err);
  protected:
    virtual int parseOptions(){ return 0;};
};

