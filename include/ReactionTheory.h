#pragma once
#include <yaml-cpp/yaml.h>
#include "TermData.h"

/**
  @class ReactionTheory

  @brief A base class manages for reaction theories

  It provides an interface which must be present in the derived classes
  Each concrete instance of ReactionTheory is a singleton
  ReactionTheory is responsible for:
  1. Getting and checking sanity of reaction parameters.
  2. Calculating theory predictions for given dataset

  @author A.Sapronov <sapronov@ifh.de>, I.Novikov<ivan.novikov@desy.de>

  @version 0.2, but nobody keeps track of versions anymore
  @date 2016/01/21
  */

//why is this not in namespace xfitter? --Ivan

class ReactionTheory{
public:
  ReactionTheory() {};
  virtual ~ReactionTheory() {};
public:
  using super=ReactionTheory;
  virtual string getReactionName()const=0; ///< Returns expected reaction name. Normally generated automatically by AddReaction.py
  virtual void atStart();            //called once after everything else is initialized
  virtual void atIteration();        //called in the beginning of each chi2 evaluation
  virtual void initTerm  (TermData*);//called once for each term, after atStart()
  virtual void reinitTerm(TermData*);//called when some parameters for this term have changed and need to be re-read
  virtual void freeTerm  (TermData*);//called for each term just before the ReactionTheory is destroyed. For cleanup
  //The following 2 methods are TEMPORARY, poorly defined and probably will be replaced
  virtual void atFCN3();
  virtual void atMakeErrorBands(int i);
  //! Main function to compute predictions for given term. Return results by filling val and errors
  virtual void compute(TermData*,valarray<double>&val,map<string,valarray<double> >&errors)=0;
protected:
  int _ncpu; // number of parallel threads
  bool _flagComputeAtIteration = false;
};
