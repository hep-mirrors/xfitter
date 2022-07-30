#pragma once

//?#include "ReactionBaseDISNC.h"
#include "ReactionTheory.h"
//#include <IntegrateDIS.h>



/*
 *   @class' ReactionFONLL_DISCC
 *
 *  @brief A wrapper class for FONLL_DISCC reaction
 *
 *  Based on the ReactionTheory class. Reads options produces 3d cross section.
 *
 *  @version 0.1
 *  @date 2017-11-29
 */

class ReactionNC_SIA : public ReactionTheory //ReactionBaseDISNC
{
 public:
  ReactionNC_SIA() {};
  //~ReactionFONLL_DISCC() {};
  //~ReactionFONLL_DISCC(const ReactionFONLL_DISCC &) {};
  //ReactionFONLL_DISCC & operator = (const ReactionAFONLL_DISCC &r) { return *(new ReactionFONLL_DISCC(r)); };

  virtual string getReactionName() const { return "NC_SIA"; };
  virtual void atStart() override final;
  virtual void atIteration() override final;
  virtual void initTerm(TermData *td)override final;
  //virtual void compute(TermData *, valarray<double> &val, map<string, valarray<double>> &errors) override final;
  virtual void compute(TermData *, valarray<double> &val, map<string, valarray<double>> &err) override;
  //virtual void initTerm(TermData*)override final;
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
  //virtual void OBSF(TermData *td, valarray<double> &val, map<string, valarray<double>> &err);
  vector<unsigned> _dsIDs; //!< list of termIDs managed by the reaction.
  map<unsigned, TermData*> _tdDS; //! store term data for later access. 


  //virtual valarray<double> FL(TermData*td) override;
  //virtual valarray<double> xF3(TermData*td) override;
  //vector<unsigned> _dsIDs;
  //map<unsigned, TermData*> _tdDS;
  //virtual const valarray<double> *GetBinValues(TermData *td, const string &binName);
 private:
  map <unsigned,valarray<double>> _obs;
  map<unsigned, dataObs> _dataObs;   //!< flavour (incl, Nor, c, b, light)
  map<unsigned, int> _npoints;         //!< Number of points in a dataset.
  map<unsigned, double> _polarisation; //!< longitudinal polarisation
  map<unsigned, double> _charge;       //!< lepton beam charge


//  map <unsigned,valarray<double>> _flfonll;
//  map <unsigned,valarray<double>> _f3fonll;
 // map <unsigned,TermData*> _dsIDs;
 protected:
  const dataObs GetDataObs(unsigned termID) { return _dataObs[termID]; }
  const double GetCharge(unsigned termID) { return _charge[termID]; }
  const int GetNpoint(unsigned termID) { return _npoints[termID]; }
  const double GetPolarisation(unsigned termID) { return _polarisation[termID]; }
  TermData* GetTermData(unsigned termID) { return _tdDS[termID];}
  //virtual const valarray<double> *GetBinValues(TermData *td, const string &binName); //! interface for integerated sigma


  bool _non_apfel_evol{false};
};
//  struct ReactionData
//  {
//    int _npoints;
//    IntegrateDIS *_integrated = nullptr;
//  };
//  const valarray<double> *GetBinValues(TermData *td, const string &binName)
//  {
//    IntegrateDIS *iDIS = ((BaseSIA::ReactionData *)td->reactionData)->_integrated;
//    if (iDIS == nullptr)
//      return &td->getBinColumn(binName);
//    if (binName == "Q2")
//      return iDIS->getBinValuesQ2();
//    else if (binName == "x")
//      return iDIS->getBinValuesX();
//    return &td->getBinColumn(binName);

//  };
//}

