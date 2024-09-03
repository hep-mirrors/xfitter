#include "ReactionTheory.h"
#include "ReactionBaseDISNC.h"
#include <IntegrateDIS.h>
#include <hoppet_v1.h>

/**
  @class' ReactionHOPPET_DISNC

  @brief HOPPET DISNC reaction

  @version 0.1
  @date 2024-08-13
  */

class ReactionHOPPET_DISNC : public ReactionBaseDISNC
{
public:
   ReactionHOPPET_DISNC(){};

public:
    virtual string getReactionName() const override { return "HOPPET_DISNC"; };

  private:
	 void F2(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal);
	 void FL(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal);
	 void xF3(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal);
    void atIteration();
};

