#pragma once

#include "ReactionTheory.h"



/*
 *   @class' ReactionNC_SIA
 *
 *  @brief A wrapper class for NC_SIA reaction
 *
 *  Based on the ReactionTheory class. Reads options produces 3d cross section.
 *
 *  @version 0.1
 *  @date 2022-5-10
 */

class ReactionNC_SIA : public ReactionTheory 
{
 public:
  ReactionNC_SIA() {};

  virtual string getReactionName() const override { return "NC_SIA"; };
  virtual void atStart() override final;
  virtual void atIteration() override final;
  virtual void initTerm(TermData *td) override final;
  virtual void compute(TermData *, valarray<double> &val, map<string, valarray<double>> &err) override;
 protected:
  enum class dataObs
  {
    incl,
    inclNor,
    inclNF4,
    inclNorNF4,
    inclCH,
    inclBT,
    inclLi
  }; //!< Define final state.
 protected:
  virtual valarray<double> OBS(TermData*td);
  vector<unsigned> _dsIDs; //!< list of termIDs managed by the reaction.
  map<unsigned, TermData*> _tdDS; //! store term data for later access. 

 private:
  map <unsigned,valarray<double>> _obs;
  map<unsigned, dataObs> _dataObs;   //!< flavour (incl, Nor, c, b, light)
  map<unsigned, int> _npoints;         //!< Number of points in a dataset.
  map<unsigned, double> _polarisation; //!< longitudinal polarisation
  map<unsigned, double> _charge;       //!< lepton beam charge


 protected:
  const dataObs GetDataObs(unsigned termID) { return _dataObs[termID]; }
  const double GetCharge(unsigned termID) { return _charge[termID]; }
  const int GetNpoint(unsigned termID) { return _npoints[termID]; }
  const double GetPolarisation(unsigned termID) { return _polarisation[termID]; }
  TermData* GetTermData(unsigned termID) { return _tdDS[termID];}


  bool _non_apfel_evol{false};
};

