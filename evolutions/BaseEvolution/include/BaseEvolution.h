#pragma once

#include <string>
#include <map>
#include <functional>

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
    BaseEvolution(const std::string& name, std::function<std::map<int,double>(double const& x)> const& inPDFs): _name(name), _inPDFs(inPDFs) { }

    /// Explicitly set PDF decomposition
    void SetPdfDecomposition(  std::function<std::map<int,double>(double const& x)> const& inPDFs) { _inPDFs = inPDFs;}
    
    /**
     * @brief Function to be called at the begining to initialise the
     * evolution code.
     */
    virtual void initAtStart() = 0;

    /**
     * @brief Function to be call at each iteration to update the
     * relevant evolution parameters.
     */
    virtual void initAtIteration() = 0;

    /**
     * @name Setters
     */
    ///@{
    /**
     * @brief Function to set the PdfDecomposition object to be used
     * as initial scale distributions.
     */
    void SetInitialPDFs(std::function<std::map<int,double>(double const& x)> const& inPDFs) { _inPDFs = inPDFs; }
    ///@}

    /**
     * @name Getters
     */
    ///@{
    /**
     * @brief Function that returns a std::function that in turns
     * returns a map<int, double> as a function of x and Q.
     * @return map<int, double>-valued function of x and Q.
     */
    virtual std::function<std::map<int,double>(double const& x, double const& Q)> xfxQMap() = 0;

    /**
     * @brief Function that returns a std::function that in turns
     * returns a double as a function of the pdf index i, x and Q.
     * @return double-valued function of i, x and Q.
     */
    virtual std::function<double(int const& i, double const& x, double const& Q)> xfxQDouble() = 0;

    /**
     * @brief Function that returns a std::function that in turns
     * returns a void as a function of the pdf index x, Q, and pdfs,
     * where pdfs is the array of PDFs.
     * @return void-valued function of x, Q and pdfs.
     */
    virtual std::function<void(double const& x, double const& Q, double* pdfs)> xfxQArray() = 0;

    /**
     * @brief Purely virtual function to get the strong coupling
     * @param Q: value of of the renormalisation scale in GeV
     * @return the strong coupling at Q.
     */
    virtual std::function<double(double const& Q)> AlphaQCD() = 0;

    const std::string getName() const { return _name; }
    
    ///@}
   
    
  protected:
    const std::string                                    _name;
    std::function<std::map<int,double>(double const& x)> _inPDFs;
  };

  /// For dynamic loader
  typedef BaseEvolution* create_evolution();
}

