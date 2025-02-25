#pragma once
#include "ReactionTheory.h"
#include "IntegrateDIS.h"
#include <spline.h>

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

   enum class dataType
   {
      signonred,
      sigred,
      f2,
      fl,
      f3,
      sigred_nof3
   }; //!< Define compute output.

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
  //const dataType GetDataType(unsigned termID) { return _dataType[termID]; }

  // higher twist
  tk::spline _spline_f2;
  tk::spline _spline_ft;
  map<unsigned, bool> _flag_ht;
  map<unsigned, std::vector<const double* > > _ht_x;
  map<unsigned, std::vector<const double* > > _ht_2;
  map<unsigned, std::vector<const double* > > _ht_t;
  map<unsigned, const double* > _ht_alpha_2;
  map<unsigned, const double* > _ht_alpha_t;
  void update_ht(size_t dataSetID);
  void apply_ht(const int dataSetID, const double q2, const double x, double& f2, double& fl);

  // nuclear corrections
  map<unsigned, int> _nucl_kint;
  map<unsigned, int> _nucl_ftyp;
  map<unsigned, int> _nucl_kord;
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
  enum class dataType
  {
    signonred,
    sigred,
    f2,
    fl,
    f3,
    sigred_nof3
  }; //!< Define compute output.
  struct ReactionData
  {
    int _npoints;                        //!< Number of points in a dataset.
    double _polarisation = 0.;           //!< longitudinal polarisation
    double _charge = 0.;                 //!< lepton beam charge
    bool _isReduced = false;             //!< reduced cross section
    bool _isBeamNu = false;              //!< neutrino beam
    dataFlav _dataFlav = dataFlav::incl; //!< flavour (incl, c, b)
    map<unsigned, dataType> _dataType;   //!< cross section (reduced, F2, FL)
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
    int _nuke_ftyp = 0; // nuclear correction: 1 Carbon (A=12, Z=6), 2 Iron (A=56, Z=26), 3 Lead (A=207, Z=82), 4 Emulsion (mixture of nuclei) 
    int _nuke_kint = 0; // nuclear interaction: 1 charged lepton, 2 neutrino (-2 antineutrino) CC, 3 neutrino (-3 antineutrino) NC, 4 neutrino CC charm production (-4 antineutrino), 5 neutrino CC not-charm (-5 antineutrino) [4, 5 not implemented]
    int _nuke_kord = 0; // order: 1 LO, 2 NLO, 3 NNLO [not implemented]
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
