#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionHathorSingleTop

  @brief A wrapper class for HathorSingleTop reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2019-06-04
// Authors: Laia Parets Peris <laia.parets.peris@desy.de>, Katerina Lipka <katerina.lipka@desy.de>
// transition from pole to MSbar scheme by S. Moch (private communication)
  */

class HathorSgTopT;
class HathorPdfxFitter;
class ReactionHathorSingleTop : public ReactionTheory
{
public:
  ReactionHathorSingleTop();
  ~ReactionHathorSingleTop();

  virtual string getReactionName() const { return  "HathorSingleTop" ;};
  virtual void initTerm(TermData *td) override final;
  virtual void atStart();
  virtual void compute(TermData *td, valarray<double> &val, map<string, valarray<double> > &err);
protected:
  virtual int parseOptions(){ return 0;};

  // this is map of key = dataset, value = pointer to Hathor instances,
  // one instance per one dataset
  std::map<int, HathorSgTopT*> _hathorArray;

  HathorPdfxFitter* _pdf;
  int* _rndStore;
  map<unsigned, int> _scheme;
  map<unsigned, double> _mtop;
  map<unsigned, double> _mr;
  map<unsigned, double> _mf;
  // store term data for later access
  map<unsigned, TermData*> _tdDS;
};

