
#pragma once

#include "ReactionTheory.h"

/**
  @class' ReactionDYTurbo

  @brief A wrapper class for DYTurbo reaction

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2019-10-07
  */

class ReactionDYTurbo : public ReactionTheory
{
public:
  ReactionDYTurbo(){};

public:
  virtual string getReactionName() const { return  "DYTurbo" ;};
  virtual void atStart();
  //  virtual void initTerm(TermData* td)override final;
  //  virtual void freeTerm(TermData* td)override final;
  virtual void compute(TermData*,valarray<double>&,map<string,valarray<double> >&)override final;
protected:
  virtual int parseOptions(){ return 0;};
};

