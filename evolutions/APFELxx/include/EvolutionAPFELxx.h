#pragma once

#include "BaseEvolution.h"

/*
 *  @class' ReactionFONLL_DISNC
 *
 *  @brief A wrapper class for FONLL_DISNC reaction 
 *
 *  Based on the ReactionTheory class. Reads options produces 3d cross section.
 *
 *  @version 0.1
 *  @date 2017-11-29
 */

namespace xfitter
{
  class EvolutionAPFELxx: public BaseEvolution
  {
  public:
    EvolutionAPFELxx(): BaseEvolution("APFELxx") {}

    /**
     * @brief Function that initialise the evolution in APFEL++.
     */
    void initAtStart() const;

  private:

  };
}
