#pragma once

#include <string>
#include <map>

#include "BasePdfDecomposition.h"

namespace xfitter
{
  /**
     @class BaseEvolution

     @brief Base class for the evolving quantities

     Contains methods to compute the evolution of PDFs, alpha_s, and
     other possible evolving quantites.

     @version 0.1
     @date 2018-07-11
  */

  class BaseEvolution
  {
  public:
    /**
     * @brief The BaseEvolution default constructor.
     * @param name: the name assignet to the instance
     */
    BaseEvolution(const std::string& name): _name(name), _PdfDecomp(nullptr) { }
    
    /**
     * @brief The BaseEvolution default constructor.
     * @param name: the name assignet to the instance
     * @param PdfDecomp: the BasePdfDecomposition object that contains distributions at the initial scale
     */
    BaseEvolution(const std::string& name, const BasePdfDecomposition& PdfDecomp): _name(name), _PdfDecomp(PdfDecomp) { }

    /**
     * @brief Function to be called at the begining to initialise the
     * evolution code.
     */
    virtual void initAtStart() const = 0;

    /**
     * @brief Function to be call at each iteration to update the
     * relevant evolution parameters.
     */ 
    virtual void initAtIteration() const = 0;

    /**
     * @name Setters
     */
    ///@{
    /**
     * @brief Function to set the PdfDecomposition object to be used
     * as initial scale distributions.
     */
    void SetPdfDecompositionObject(BasePdfDecomposition const& PdfDecomp) { _PdfDecomp = PdfDecomp; }
    ///@}

    /**
     * @name Getters
     */
    ///@{
    /**
     * @brief Function to get evolved distributions
     * @param x: value of Bjorken n
     * @param Q: value of of the factorisation scale in GeV
     * @return a map of int to double containing the relevant distributions at x and Q.
     */
    virtual std::map<int, double> const& GetDistributions(double const& x, double const& Q) const = 0;

    /**
     * @brief Purely virtual function to get the strong coupling
     * @param Q: value of of the renormalisation scale in GeV
     * @return the strong coupling at Q.
     */
    virtual double const& GetAlphaQCD(double const& Q) const = 0;
    ///@}

  private:
    /// Name of the evolution object
    const std::string _name;

    /// Input PDF set
    BasePdfDecomposition _PdfDecomp;
  };
}
