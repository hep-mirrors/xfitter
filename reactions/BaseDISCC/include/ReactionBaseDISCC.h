
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
#define BASE_PARS (int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)

class ReactionBaseDISCC : public ReactionTheory
{
  public:
    ReactionBaseDISCC(){};

//    ~ReactionBaseDISCC(){};
//    ~ReactionBaseDISCC(const ReactionBaseDISCC &){};
//    ReactionBaseDISCC & operator =(const ReactionABaseDISCC &r){return *(new ReactionBaseDISCC(r));};

  public:
    virtual string getReactionName() const { return  "BaseDISCC" ;};
    int initAtStart(const string &); 
    virtual void setDatasetParamters( int dataSetID, map<string,string> pars, map<string,double> parsDataset) override ;
    virtual void initAtIteration() override;
    
    virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err);
  protected:
    virtual int parseOptions(){ return 0;};

    virtual void F2 BASE_PARS;
    virtual void FL BASE_PARS;
    virtual void xF3 BASE_PARS;

 private:
    map <int, int>     _npoints ;           //!< Number of points in a dataset.
    map <int, double>  _polarisation ;      //!< longitudinal polarisation
    map <int, double>  _charge;             //!< lepton beam charge
    map <int, int>    _isReduced;          //!< reduced cross section
  protected:
    const int GetNpoint(int dataSetID) {return _npoints[dataSetID];}
    const double GetPolarisation (int dataSetID) {return _polarisation[dataSetID];}
    const double GetCharge(int dataSetID) {return _charge[dataSetID]; }
    const int IsReduced(int dataSetID){ return _isReduced[dataSetID] > 0; }
    
    // Another decomposition:
    virtual void GetF2u( int dataSetID, valarray<double>& f2u);
    virtual void GetFLu( int dataSetID, valarray<double>& flu);
    virtual void GetxF3u( int dataSetID, valarray<double>& xf3u );
    virtual void GetF2d( int dataSetID, valarray<double>& f2d);
    virtual void GetFLd( int dataSetID, valarray<double>& fld);
    virtual void GetxF3d( int dataSetID, valarray<double>& xf3d );
    
  private:
    // Some buffering mechanism to avoid double calls
    map <int,valarray<double> > _f2u; //!< F2 for u-type quarks
    map <int,valarray<double> > _f2d; //!< F2 for d-type quarks
    map <int,valarray<double> > _flu; //!< FL for u-type quarks
    map <int,valarray<double> > _fld; //!< FL for d-type quarks
    map <int,valarray<double> > _xf3u; 
    map <int,valarray<double> > _xf3d;

 protected:
    double _MW;
    double _Gf;
    double _convfac;
};

