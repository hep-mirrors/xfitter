
#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionBaseDISCC

  @brief A wrapper class for BaseDISCC reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2017-10-05
  */

// Define standard parameters used by SF and x-sections:
// What?? --Ivan
#define BASE_PARS (int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)

class IntegrateDIS;

class ReactionBaseDISCC : public ReactionTheory
{
  public:
    ReactionBaseDISCC(){};
  public:
    virtual string getReactionName() const { return  "BaseDISCC" ;};
    virtual void atStart()override final;
    virtual void atIteration()override final;
    virtual void initTerm(TermData*)override final;
    virtual void compute(TermData*,valarray<double>&val,map<string,valarray<double> >&errors)override final;
  protected:
    virtual void F2 BASE_PARS;
    virtual void FL BASE_PARS;
    virtual void xF3 BASE_PARS;

 private:
 //All this is deep WIP!! --Ivan
  protected:
    const int GetNpoint(int dataSetID) {return _npoints[dataSetID];}
    const double GetPolarisation (int dataSetID) {return _polarisation[dataSetID];}
    const double GetCharge(int dataSetID) {return _charge[dataSetID]; }
    const int IsReduced(int dataSetID){ return _isReduced[dataSetID] > 0; }
    const dataFlav GetDataFlav(int dataSetID) {return _dataFlav[dataSetID]; }

    // Another decomposition:
    virtual void GetF2u( int dataSetID, valarray<double>& f2u);
    virtual void GetFLu( int dataSetID, valarray<double>& flu);
    virtual void GetxF3u( int dataSetID, valarray<double>& xf3u );
    virtual void GetF2d( int dataSetID, valarray<double>& f2d);
    virtual void GetFLd( int dataSetID, valarray<double>& fld);
    virtual void GetxF3d( int dataSetID, valarray<double>& xf3d );

  protected:
    virtual valarray<double> *GetBinValues(int idDS, const string& binName);

 protected:
    double _MW;
    double _Gf;
    double _convfac;
};

