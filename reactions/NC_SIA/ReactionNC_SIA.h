#pragma once

#include "ReactionBaseDISCC.h"

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

class ReactionNC_SIA : public ReactionBaseDISCC
{
 public:
  ReactionNC_SIA() {};
  //~ReactionFONLL_DISCC() {};
  //~ReactionFONLL_DISCC(const ReactionFONLL_DISCC &) {};
  //ReactionFONLL_DISCC & operator = (const ReactionAFONLL_DISCC &r) { return *(new ReactionFONLL_DISCC(r)); };

  virtual string getReactionName() const { return "NC_SIA"; };
  virtual void atStart() override final;
  virtual void atIteration() override final;
  virtual void initTerm(TermData*) override final;
  //virtual void compute(TermData *, valarray<double> &val, map<string, valarray<double>> &errors) override final;

  //virtual void initTerm(TermData*)override final;
 
 protected:
  virtual valarray<double> OBS(TermData*td);
  //virtual valarray<double> FL(TermData*td) override;
  //virtual valarray<double> xF3(TermData*td) override;
  //vector<unsigned> _dsIDs;
  //map<unsigned, TermData*> _tdDS;
  //virtual const valarray<double> *GetBinValues(TermData *td, const string &binName);
// private:
  map <unsigned,valarray<double>> _obs;
//  map <unsigned,valarray<double>> _flfonll;
//  map <unsigned,valarray<double>> _f3fonll;
  map <unsigned,TermData*> _dsIDs;
};
//namespace BaseSIA
//{
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

