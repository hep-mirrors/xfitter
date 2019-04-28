#pragma once
#include "ReactionTheory.h"

/**
  @class' ReactionBaseDISCC

  @brief A wrapper class for BaseDISCC reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-10-05
  */
class ReactionBaseDISCC : public ReactionTheory
{
  public:
    ReactionBaseDISCC(){};
  public:
    virtual string getReactionName()const{return"BaseDISCC";};
    virtual void atStart()override final;
    virtual void atIteration()override final;
    virtual void initTerm(TermData*)override final;
    virtual void compute(TermData*,valarray<double>&val,map<string,valarray<double> >&errors)override final;
  private:
    //these 5 were used in derived classes
    //TODO: reimplement derived classes
/*
    const int GetNpoint(int dataSetID) {return _npoints[dataSetID];}
    const double GetPolarisation (int dataSetID) {return _polarisation[dataSetID];}
    const double GetCharge(int dataSetID) {return _charge[dataSetID]; }
    const int IsReduced(int dataSetID){ return _isReduced[dataSetID] > 0; }
    const dataFlav GetDataFlav(int dataSetID) {return _dataFlav[dataSetID]; }
*/
  protected:
    double _Gf;
    double _convfac;
};

