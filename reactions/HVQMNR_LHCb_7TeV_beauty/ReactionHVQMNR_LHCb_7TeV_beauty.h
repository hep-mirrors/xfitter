
#pragma once

#include "ReactionBaseHVQMNR.h"

/**
  @class' ReactionHVQMNR_LHCb_7TeV_beauty

  @brief A wrapper class for HVQMNR_LHCb_7TeV_beauty reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.
  *
  Derived from ReactionBaseHVQMNR where basic stuff for HVQMNR calculation is implemented.
  This class implements calculation for LHCb beauty measurement at 7 TeV
  [JHEP 1308 (2013) 117] [arXiv:1306.3663]

  @version 0.1
  @date 2017-01-02
  */

class ReactionHVQMNR_LHCb_7TeV_beauty : public ReactionBaseHVQMNR
{
public:
  ReactionHVQMNR_LHCb_7TeV_beauty(){};

public:
  virtual string getReactionName() const override { return  "HVQMNR_LHCb_7TeV_beauty" ;};
  virtual void initTerm(TermData *td) override final;
  virtual void atIteration() override final;        //called in the beginning of each chi2 evaluation
  virtual void compute(TermData *td, valarray<double> &val, map<string, valarray<double>> &errors) override final;
protected:
  virtual int parseOptions(){ return 0;};
};

