#pragma once
#include "ReactionTheory.h"
#include "IntegrateDIS.h"

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
  virtual string getReactionName() const override { return "BaseDISCC"; };
  virtual void atStart() override;
  virtual void initTerm(TermData *) override;
  virtual void compute(TermData *, valarray<double> &val, map<string, valarray<double>> &errors) override final;

protected:

  virtual valarray<double> FL(TermData *td);
  virtual valarray<double> F2(TermData *td);
  virtual valarray<double> xF3(TermData *td);

  double _Gf;
  double _convfac;

  vector<unsigned> _dsIDs; //!< list of termIDs managed by the reaction.
  map<unsigned, TermData*> _tdDS; //! store term data for later access.

  // for integrated cross sections
  // method is based on legacy subroutine GetIntegratedDisXsection
  // is the map with pointers IntegrateDIS* outdated (seems to be in ReactionData)? --SZ
  //map<unsigned, IntegrateDIS *> _integrated;
  virtual const valarray<double> *GetBinValues(TermData *td, const string &binName); //! interface for integerated sigma
};

/// Helper classes and functions:
namespace BaseDISCC
{
  enum class dataFlav
  {
    incl,
    c,
    b
  }; //!< Define final state.
  struct ReactionData
  {
    int _npoints;                        //!< Number of points in a dataset.
    double _polarisation = 0.;           //!< longitudinal polarisation
    double _charge = 0.;                 //!< lepton beam charge
    bool _isReduced = false;             //!< reduced cross section
    dataFlav _dataFlav = dataFlav::incl; //!< flavour (incl, c, b)
    // for integrated cross sections
    // method is based on legacy subroutine GetIntegratedDisXsection
    IntegrateDIS *_integrated = nullptr;
    // Some buffering mechanism to avoid double calls
    valarray<double> _f2u; //!< F2 for u-type quarks
    valarray<double> _f2d; //!< F2 for d-type quarks
    valarray<double> _flu; //!< FL for u-type quarks
    valarray<double> _fld; //!< FL for d-type quarks
    valarray<double> _xf3u;
    valarray<double> _xf3d;
    const double *Mw; //parameter of W mass
    int ipdfSet; /// PDF set used in the QCDNUM evolution
  };

  /// Helper function to get bin values, including integrated sigma:
  const valarray<double> *GetBinValues(TermData *td, const string &binName)
  {
    IntegrateDIS *iDIS = ((BaseDISCC::ReactionData *)td->reactionData)->_integrated;
    if (iDIS == nullptr)
      return &td->getBinColumn(binName);
    if (binName == "Q2")
      return iDIS->getBinValuesQ2();
    else if (binName == "x")
      return iDIS->getBinValuesX();
    else if (binName == "y")
      return iDIS->getBinValuesY();
    return &td->getBinColumn(binName);
  };
} // namespace BaseDISCC
