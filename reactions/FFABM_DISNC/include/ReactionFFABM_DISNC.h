
#pragma once

#include "ReactionBaseDISNC.h"

/**
  @class' ReactionFFABM_DISNC

  @brief A wrapper class for FFABM_DISNC reaction 

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-09-29
  */

class ReactionFFABM_DISNC : public ReactionBaseDISNC
{
  private:
    typedef ReactionBaseDISNC Super;
  public:
    ReactionFFABM_DISNC(){};
  public:
    virtual string getReactionName() const { return  "FFABM_DISNC" ;};
    int initAtStart(const string &);
    virtual void setDatasetParamters( int dataSetID, map<string,string> pars, map<string,double> parsDataset) override ;
    virtual void initAtIteration() override;

  protected:
    virtual void F2 BASE_PARS override;
    virtual void FL BASE_PARS override;
    virtual void xF3 BASE_PARS override;

  private:
    map <int,valarray<double> > _f2abm;
    map <int,valarray<double> > _flabm;
    map <int,valarray<double> > _f3abm;

    // parameters initialised at iteration
    double _mc;
    double _mb;
    double _mz;
    double _asmz;
    double _sin2thw;
    double _cos2thw;

    void calcF2FL(int dataSetID);
};

