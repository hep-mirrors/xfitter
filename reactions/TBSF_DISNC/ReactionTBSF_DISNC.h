#pragma once

#include "ReactionBaseDISNC.h"
#include "ReactionBaseTBSF.h"
#include "ForkPool.h"

class ReactionTBSF_DISNC : public ReactionBaseTBSF<ReactionBaseDISNC, abm::SFproc::nc>
{
public:
  virtual string getReactionName() const override { return "TBSF_DISNC"; };

protected:
  virtual void F2 BASE_PARS override final { val = _f2abm[td->id]; };
  virtual void FL BASE_PARS override final { val = _flabm[td->id]; };
  virtual void xF3 BASE_PARS override final { val = _f3abm[td->id]; };
};
