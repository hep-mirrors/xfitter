
#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionTensorPomeron

  @brief A wrapper class for TensorPomeron reaction 

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-04-14
  */

class ReactionTensorPomeron : public ReactionTheory
{
  public:
    ReactionTensorPomeron(){};

//    ~ReactionTensorPomeron(){};
//    ~ReactionTensorPomeron(const ReactionTensorPomeron &){};
//    ReactionTensorPomeron & operator =(const ReactionATensorPomeron &r){return *(new ReactionTensorPomeron(r));};

  public:
    virtual string getReactionName() const { return  "TensorPomeron" ;};
    int initAtStart(const string &); 
    virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err);
  protected:
    virtual int parseOptions(){ return 0;};
};

