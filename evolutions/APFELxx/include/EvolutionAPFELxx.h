#pragma once

#include "BaseEvolution.h"

#include <apfel/grid.h>
#include <apfel/dglap.h>
#include <apfel/dglapbuilder.h>
#include <apfel/tabulateobject.h>

#include <vector>
#include <memory>
#include <functional>

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
    /// Empty constructor (needed for the dynamic loading)
    EvolutionAPFELxx():  BaseEvolution{"APFELxx",nullptr} {}

    /// Constructor wit PDF decomposition
    EvolutionAPFELxx(std::function<std::map<int,double>(double const& x)> const& inPDFs): BaseEvolution{"APFELxx", inPDFs} {}

    /**
     * @brief Function that initialises the evolution in APFEL++.
     */
    void initAtStart();

    /**
     * @brief Function that updates the relevant parameters of the
     * evolution at each iteration of the fitting procedure.
     */
    void initAtIteration();

    /**
     * @name Getters
     */
    ///@{
    /**
     * @brief Function that returns a std::function that in turns
     * returns a map<int, double> as a function of x and Q.
     * @return map<int, double>-valued function of x and Q.
     */
    std::function<std::map<int,double>(double const& x, double const& Q)> xfxQMap();

    /**
     * @brief Function that returns a std::function that in turns
     * returns a double as a function of the pdf index i, x and Q.
     * @return double-valued function of i, x and Q.
     */
    std::function<double(int const& i, double const& x, double const& Q)> xfxQDouble();

    /**
     * @brief Function that returns a std::function that in turns
     * returns a void as a function of the pdf index x, Q, and pdfs,
     * where pdfs is the array of PDFs.
     * @return void-valued function of x, Q and pdfs.
     */
    std::function<void(double const& x, double const& Q, double* pdfs)> xfxQArray();

    /**
     * @brief Function that returns a std::function that in turns
     * returns double as a function of  Q.
     * @return double-valued function of Q.
     */
    std::function<double(double const& Q)> AlphaQCD() { return _AlphaQCD; };
    ///@}

  private:
    std::vector<double>                                                     _Masses;
    std::vector<double>                                                     _Thresholds;
    std::unique_ptr<const apfel::Grid>                                      _Grid;
    std::map<int,apfel::DglapObjects>                                       _DglapObj;
    std::unique_ptr<apfel::TabulateObject<apfel::Set<apfel::Distribution>>> _TabulatedPDFs;
    std::function<double(double const& Q)>                                  _AlphaQCD;
  };
}
