#pragma once

#include "ReactionTheory.h"
#include "ReactionBaseDISNC.h"
#include <IntegrateDIS.h>

/**
  @class' ReactionBaseDISNC

  @brief A wrapper class for BaseDISNC reaction

  Based on the ReactionTheory class.

  @version 0.1
  @date 2017-04-08
  */

// Define standard parameters used by SF and x-sections:
#define BASE_PARS (TermData * td, valarray<double> & val, map<string, valarray<double>> & err)

class Reaction_DISNC_Hoppet : public ReactionBaseDISNC
{
public:
   Reaction_DISNC_Hoppet(){};

public:
   virtual string getReactionName() const override { return "HoppetDISNC"; };
   virtual void atStart() override;
   virtual void initTerm(TermData *td) override;

   //!< Initialize all EWK couplings here:
   virtual void atIteration() override;
   virtual void compute(TermData *, valarray<double> &val, map<string, valarray<double>> &err) override;

protected:
   /*
       A few methods specific for DIS NC process.
    */

   virtual void F2gamma BASE_PARS override;
   virtual void F2gammaZ BASE_PARS override;
   virtual void F2Z BASE_PARS override ;


   virtual void FLgamma BASE_PARS override ;
   virtual void FLgammaZ BASE_PARS override;
   virtual void FLZ BASE_PARS override;

   virtual void xF3gammaZ BASE_PARS override;
   virtual void xF3Z BASE_PARS override ;
};
