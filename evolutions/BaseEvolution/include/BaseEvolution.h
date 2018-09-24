#pragma once

#include <string>
#include <map>
#include <functional>
#include <yaml-cpp/yaml.h>

namespace xfitter
{
  /**
     @class BaseEvolution

     @brief Base class for the evolved PDFs and other evolved quantities

     Contains methods to compute the evolution of PDFs, alpha_s, and
     other possible evolving quantites.

     @version 0.2
     @date 2018-09-29
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
    virtual void initAtStart(){};

    /**
     * @brief Function to be called at each iteration
     */
    virtual void initAtIteration(){};

    /**
     * @brief This function should be called when at least one parameter in the YAML node of given evolution changes
     */
    virtual void initAtParameterChange(){};

    /**
     * @name Getters
     */
    ///@{
    /**
     * @brief Function that returns a std::function that in turn
     * returns a map<int, double> as a function of x and Q.
     * @return map<int, double>-valued function of x and Q.
     */
    virtual std::function<std::map<int,double>(double const& x, double const& Q)> xfxQMap() = 0;

    /**
     * @brief Function that returns a std::function that in turn
     * returns a double as a function of the pdf index i, x and Q.
     * @return double-valued function of i, x and Q.
     */
    virtual std::function<double(int const& i, double const& x, double const& Q)> xfxQDouble() = 0;

    /**
     * @brief Function that returns a std::function that in turn
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

    /// Get generic property of the evolution
    virtual std::string getPropertyS(std::string const& propertyName ) const { return "" ; }

    /// Get generic property of the evolution
    virtual int  getPropertyI(std::string const& propertyName ) const { return 0; }

    /// Get generic property of the evolution
    virtual double  getPropertyD(std::string const& propertyName ) const { return 0.; }
    
    /// Get class name, can be used to verify that the correct concrete class is being used
    virtual const char*getClassName()const=0;
    ///@}
   
    
  protected:
    std::function<std::map<int,double>(double const& x)> _inPDFs;
  };

  /// For dynamic loader
  typedef BaseEvolution*create_evolution(const char*name);
}

