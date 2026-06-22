
#pragma once

#include "ReactionBaseDISCC.h"
#include "ReactionBaseFFABM.h"
#include "ForkPool.h"

class ReactionFFABM_DISNC_CC : public ReactionBaseFFABM, public ReactionBaseDISCC
{
private:
  typedef ReactionBaseDISCC Super;

public:
  virtual string getReactionName() const override { return "FFABM_DISNC_CC"; };
  virtual void initTerm(TermData *td) { Super::initTerm(td); ReactionBaseFFABM::initTerm(td); };
  virtual void atIteration() override { Super::atIteration(); ReactionBaseFFABM::atIteration(); };

private:
  virtual abm::SFproc GetProc() { return abm::SFproc::cc; };
  virtual void compute(TermData *td, valarray<double> &val, map<string, valarray<double>> &errors) {return Super::compute(td, val, errors); };
  virtual const valarray<double> *GetBinValues(TermData *td, const string &binName) { return Super::GetBinValues(td, binName); }
  virtual const double GetPolarisation(unsigned termID) { return Super::GetPolarisation(termID); };
  virtual const double GetCharge(unsigned termID)  { return Super::GetCharge(termID); };
  const ReactionBaseFFABM::dataFlav GetDataFlav(unsigned termID) { return ReactionBaseFFABM::dataFlav(Super::GetDataFlav(termID)); }

protected:
  virtual valarray<double> F2(TermData *td) override final { return _f2abm[td->id]; };
  virtual valarray<double> FL(TermData *td) override final { return _flabm[td->id]; };
  virtual valarray<double> xF3(TermData *td) override final { return _f3abm[td->id]; };
};
