
#pragma once

#include "ReactionBaseDISNC.h"
#include "ReactionBaseFFABM.h"
#include "ForkPool.h"

class ReactionFFABM_DISNC_NC : public ReactionBaseFFABM, public ReactionBaseDISNC
{
private:
  typedef ReactionBaseDISNC Super;

public:
  virtual string getReactionName() const override { return "FFABM_DISNC_NC"; };
  virtual void initTerm(TermData *td) { Super::initTerm(td); ReactionBaseFFABM::initTerm(td); };
  virtual void atIteration() override { Super::atIteration(); ReactionBaseFFABM::atIteration(); };

private:
  virtual abm::SFproc GetProc() { return abm::SFproc::nc; };
  virtual void compute(TermData *td, valarray<double> &val, map<string, valarray<double>> &errors) {return Super::compute(td, val, errors); };
  virtual const valarray<double> *GetBinValues(TermData *td, const string &binName) { return Super::GetBinValues(td, binName); }
  virtual const double GetPolarisation(unsigned termID) { return Super::GetPolarisation(termID); };
  virtual const double GetCharge(unsigned termID)  { return Super::GetCharge(termID); };
  const ReactionBaseFFABM::dataFlav GetDataFlav(unsigned termID) { return ReactionBaseFFABM::dataFlav(Super::GetDataFlav(termID)); }

protected:
  virtual void F2 BASE_PARS override final { val = _f2abm[td->id]; };
  virtual void FL BASE_PARS override final { val = _flabm[td->id]; };
  virtual void xF3 BASE_PARS override final { val = _f3abm[td->id]; };
};
