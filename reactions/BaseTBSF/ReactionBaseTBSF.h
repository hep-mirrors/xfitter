
#pragma once
#include "ReactionBaseFFABM.h"
//#include "cuba.h"

class TermData;

template <class BaseDIS, abm::SFproc Proc>
class ReactionBaseTBSF : public ReactionBaseFFABM<BaseDIS, Proc> {
public:
  //~ReactionBaseTBSF();
  //virtual void initTerm(TermData *td) override;
  //virtual void atIteration() override;

protected:
  struct DataPoint: public ReactionBaseFFABM<BaseDIS, Proc>::DataPoint {
    virtual void calc_at_q2x();
  };
};
