#pragma once

#include <string>
#include <map>
#include <vector>

namespace xfitter
{
  /**
     @class BaseEvolution

     @brief Base class for the evolved PDFs and other evolved quantities

     Contains methods to compute the evolution of PDFs, alpha_s, and
     other possible evolving quantites.

     @version 0.3
     @date 2019-04-22
  */
  class BaseEvolution
  {
  public:
    /// Unique name of instance
    const std::string _name;
    /**
     * @brief The BaseEvolution default constructor.
     * @param name: the unique name used to identify the instance
     */
    BaseEvolution(const char*name):_name(name){}

    /**
     * @brief Function to be called at the begining to initialise the
     * evolution code based on its YAML node
     *
     * This function is called only once
     */
    virtual void atStart(){};

    /**
     * @brief Function to be called at each iteration
     */
    virtual void atIteration(){};

    /**
     * @brief Function to be called at each iteration, after atIteration() has been called for all evolutions, but before theory predictions are computed
     */
    virtual void afterIteration(){};

    /**
     * @brief This function should be called when at least one parameter in the YAML node of given evolution changes
     */
    virtual void atConfigurationChange(){};

    /**
     * @brief This function should be called when writing out PBTMDs
     */
    virtual void Write_TMD(const char* name){};

    /**
     * @name Getters
     */
    ///@{
    /**
     * @brief Get PDFs for given x and Q as a map<int,double>
     * @return PDF as map<int,double>
     *Parton codes are LHAPDF convention:
     *      i  -6 -5 -4 -3 -2 -1 21 1  2  3  4  5  6
     * pdfs[i] tb bb cb sb ub db g  d  u  s  c  b  t
     */
    virtual std::map<int,double>xfxQmap(double x,double Q)=0;

    /**
     * @brief Get PDF for given x, Q, and flavor i
     * i indexes flavor (QCDNUM convention):
     *
     * i -6 -5 -4 -3 -2 -1 0  1  2  3  4  5  6
     *   tb bb cb sb ub db g  d  u  s  c  b  t
     *
     * @return PDF value
     */
    virtual double xfxQ(int i,double x,double Q)=0;

    /**
     * @brief Get PDF for given x and Q by filling the given array
     * writes into array pdfs (C++ QCDNUM convention):
     *
     *      i  0  1  2  3  4  5  6  7  8  9 10 11 12
     * pdfs[i] tb bb cb sb ub db g  d  u  s  c  b  t
     *
     * @return PDFs
     */
    virtual void xfxQarray(double x,double Q,double*pdfs)=0;

    /**
     * @brief Get the strong coupling constant at scale Q
     * @param Q: value of of the renormalisation scale in GeV
     * @return the strong coupling at Q.
     */
    virtual double getAlphaS(double Q)=0;

    /// Get generic property of the evolution
    virtual std::string getPropertyS(std::string const& propertyName ) const { return "" ; }

    /// Get generic property of the evolution
    virtual int  getPropertyI(std::string const& propertyName ) const { return 0; }

    /// Get generic property of the evolution
    virtual double  getPropertyD(std::string const& propertyName ) const { return 0.; }

    /// Get generic property of the evolution with default value and no error
    virtual double  getPropertyD(std::string const& propertyName, double defval ) const { return -1.; }

    /// Get class name, can be used to verify that the correct concrete class is being used
    virtual const char*getClassName()const=0;

    /**
     * @brief Get x points of internal grid
     * @return A sorted array of x points, or an empty array if such grid information is not available
     */
    virtual std::vector<double> getXgrid(){return std::vector<double>();}

    /**
     * @brief Get Q points of internal grid
     * @return A sorted array of Q points, or an empty array if such grid information is not available
     */
    virtual std::vector<double> getQgrid(){return std::vector<double>();}
    ///@}
  };
}

