#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionBaseDISNC

  @brief A wrapper class for BaseDISNC reaction

  Based on the ReactionTheory class.

  @version 0.1
  @date 2017-04-08
  */

// Define standard parameters used by SF and x-sections:
#define BASE_PARS (int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)

class IntegrateDIS;

class ReactionBaseDISNC : public ReactionTheory
{
 public:
    ReactionBaseDISNC(){};
 public:
    virtual string getReactionName() const { return  "BaseDISNC" ;};
    int atStart(const string &);
    virtual void setDatasetParameters( int dataSetID, map<string,string> pars, map<string,double> parsDataset) override ;

    //!< Initialize all EWK couplings here:
    virtual void initAtIteration() override;
    virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err) override ;
 protected:
    enum class dataType { signonred, sigred, f2, fl} ;  //!< Define compute output.
    enum class dataFlav { incl, c, b} ;      //!< Define final state.

    /*
       A few methods specific for DIS NC process.
    */

    virtual void F2gamma  BASE_PARS ;
    virtual void F2gammaZ BASE_PARS ;
    virtual void F2Z      BASE_PARS ;

    //! compute full F2
    virtual void F2 BASE_PARS ;

    virtual void FLgamma  BASE_PARS ;
    virtual void FLgammaZ BASE_PARS ;
    virtual void FLZ      BASE_PARS ;

    //!< compute full FL
    virtual void FL BASE_PARS ;

    virtual void xF3gammaZ BASE_PARS ;
    virtual void xF3Z      BASE_PARS ;

    //!< compute full xF3
    virtual void xF3 BASE_PARS;

    //! reduced cross section
    virtual void sred BASE_PARS;

    // Helper functions:
    void kappa(int dataSetID, valarray<double>& k) ;

 private:
    map <int, int>     _npoints ;           //!< Number of points in a dataset.
    map <int, double>  _polarisation ;      //!< longitudinal polarisation
    map <int, double>  _charge;             //!< lepton beam charge
    map <int, dataType> _dataType;          //!< cross section (reduced, F2, FL)
    map <int, dataFlav> _dataFlav;          //!< flavour (incl, c, b)

 protected:
    // some parameters which may change from iteration to iteration:
    double _alphaem;
    double _Mz;
    double _Mw;
    double _sin2thetaW;
    double _ae, _ve;
    double _au, _ad;
    double _vu, _vd;

    // conversion constant factor
    double _convfac;

 protected:
    const int GetNpoint(int dataSetID) {return _npoints[dataSetID];}
    const double GetPolarisation (int dataSetID) {return _polarisation[dataSetID];}
    const double GetCharge(int dataSetID) {return _charge[dataSetID]; }
    const dataType GetDataType(int dataSetID) {return _dataType[dataSetID]; }
    const dataFlav GetDataFlav(int dataSetID) {return _dataFlav[dataSetID]; }

    // Another decomposition:
    virtual void GetF2ud( int dataSetID, valarray<double>& f2u, valarray<double>& f2d);
    virtual void GetFLud( int dataSetID, valarray<double>& flu, valarray<double>& fld);
    virtual void GetxF3ud( int dataSetID, valarray<double>& xf3u, valarray<double>& xf3d );

 private:
    // Some buffering mechanism to avoid double calls
    map <int,valarray<double> > _f2u; //!< F2 for u-type quarks
    map <int,valarray<double> > _f2d; //!< F2 for d-type quarks
    map <int,valarray<double> > _flu; //!< FL for u-type quarks
    map <int,valarray<double> > _fld; //!< FL for d-type quarks
    map <int,valarray<double> > _xf3u;
    map <int,valarray<double> > _xf3d;

  protected:
    // for integrated cross sections
    // method is based on legacy subroutine GetIntegratedDisXsection
    map<int,IntegrateDIS*> _integrated;
    virtual valarray<double> *GetBinValues(int idDS, const string& binName);
};

