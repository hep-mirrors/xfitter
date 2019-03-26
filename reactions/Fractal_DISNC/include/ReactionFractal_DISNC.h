
#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionFractal_DISNC

  @brief A wrapper class for Fractal_DISNC reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-05-15
  */

class ReactionFractal_DISNC : public ReactionTheory
{
  public:
    ReactionFractal_DISNC(){};

//    ~ReactionFractal_DISNC(){};
//    ~ReactionFractal_DISNC(const ReactionFractal_DISNC &){};
//    ReactionFractal_DISNC & operator =(const ReactionAFractal_DISNC &r){return *(new ReactionFractal_DISNC(r));};

  public:
    virtual string getReactionName() const { return  "Fractal_DISNC" ;};
    int atStart(const string &){return 0;};
    virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err);
  protected:
//    virtual int parseOptions(){ return 0;};
};

