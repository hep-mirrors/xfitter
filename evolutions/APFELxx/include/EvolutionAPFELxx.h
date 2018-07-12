#pragma once

#include "BaseEvolution.h"

#include <apfel/grid.h>
#include <apfel/dglap.h>
#include <apfel/dglapbuilder.h>
#include <apfel/tabulateobject.h>

#include <vector>
#include <memory>

namespace xfitter
{
  /**
     @class EvolutionAPFELxx

     @brief Derived class of BaseEvolution for using APFEL++ as an evolution code.

     @version 0.1
     @date 2018-07-11
  */
  class EvolutionAPFELxx: public BaseEvolution
  {
  public:
    EvolutionAPFELxx(): BaseEvolution("APFELxx") {}

    /**
     * @brief Function that initialise the evolution in APFEL++.
     */
    void initAtStart();

    /**
     * @brief Function that updates the relevant parameters of the
     * evolution at each of the fitting procedure.
     */
    void initAtIteration();

  private:
    std::vector<double>                                                     _Masses;
    std::vector<double>                                                     _Thresholds;
    std::map<int,apfel::DglapObjects>                                       _DglapObj;
    std::unique_ptr<apfel::TabulateObject<apfel::Set<apfel::Distribution>>> _TabulatedPDFs;
  };
}
