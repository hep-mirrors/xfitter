
#pragma once

#include "ReactionBaseDISCC.h"

/**
  @class' ReactionFFABM_DISCC

  @brief A wrapper class for FFABM_DISCC reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-10-09
  */

class ReactionFFABM_DISCC : public ReactionBaseDISCC
{
  private:
    typedef ReactionBaseDISCC Super;
  public:
    ReactionFFABM_DISCC(){};
    virtual string getReactionName() const { return  "FFABM_DISCC" ;};
    int atStart(const string &);
    virtual void setDatasetParameters( int dataSetID, map<string,string> pars, map<string,double> parsDataset) override ;
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

