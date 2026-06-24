#pragma once

#include "ReactionBaseDISCC.h"
#include "ReactionBaseFFABM.h"
#include "ForkPool.h"

class ReactionFFABM_DISCC : public ReactionBaseFFABM<ReactionBaseDISCC, abm::SFproc::cc>
{
public:
  virtual string getReactionName() const override { return "FFABM_DISCC"; };

protected:
  virtual valarray<double> F2(TermData *td) override final { return _f2abm[td->id]; };
  virtual valarray<double> FL(TermData *td) override final { return _flabm[td->id]; };
  virtual valarray<double> xF3(TermData *td) override final { return _f3abm[td->id]; };
};
