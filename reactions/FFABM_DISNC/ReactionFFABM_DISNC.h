#pragma once

#include "ReactionBaseDISNC.h"
#include "ReactionBaseFFABM.h"
#include "ForkPool.h"

class ReactionFFABM_DISNC : public ReactionBaseFFABM<ReactionBaseDISNC, abm::SFproc::nc>
{
public:
  virtual string getReactionName() const override { return "FFABM_DISNC"; };

protected:
  virtual void F2 BASE_PARS override final { val = _f2abm[td->id]; };
  virtual void FL BASE_PARS override final { val = _flabm[td->id]; };
  virtual void xF3 BASE_PARS override final { val = _f3abm[td->id]; };
};
