
#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionfastNLO

  @brief A wrapper class for fastNLO reaction 

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2016-12-06
  */

class ReactionfastNLO : public ReactionTheory
{
  public:
    ReactionfastNLO(){};

//    ~ReactionfastNLO(){};
//    ~ReactionfastNLO(const ReactionfastNLO &){};
//    ReactionfastNLO & operator =(const ReactionAfastNLO &r){return *(new ReactionfastNLO(r));};

  public:
    virtual string getReactionName() const { return  "fastNLO" ;};
    int  initAtStart(const string &); 
    virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err);
  protected:
    virtual int parseOptions(){ return 0;};
};

