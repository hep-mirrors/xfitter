
#pragma once

/**
  @class' EvolutionLHAPDF

  @brief A wrapper class for LHAPDF evolution 

  @version 0.1
  @date 2018-08-20
  */

#include "BaseEvolution.h"
#include "LHAPDF/LHAPDF.h"

namespace xfitter 
{

class EvolutionLHAPDF : BaseEvolution 
{
  public:
    /// Empty constructor (needed for the dynamic loading)
    EvolutionLHAPDF():  BaseEvolution("LHAPDF",nullptr) {};

  public:
  /// Constructor setting the name
    virtual std::string getEvolutionName() const { return  "LHAPDF" ;};
  /// Global initialization
    virtual void initAtStart() override final;
  /// Init at each iteration
    virtual void initAtIteration() override final;  

  /// Return PDFs as a map <int,double> where int is PDF ID (-6, ... 6, 21)   
    virtual std::function<std::map<int,double>(double const& x, double const& Q)> xfxQMap() override final;

  /// Returns PDFs as a function of i, x, Q
    virtual std::function<double(int const& i, double const& x, double const& Q)> xfxQDouble() override final;
    
  /// Returns PDFs as double pdfs* --> double[13] from -6 to 6.  
    virtual std::function<void(double const& x, double const& Q, double* pdfs)> xfxQArray() override final;

  /// Returns alphaS
    virtual std::function<double(double const& Q)> AlphaQCD() override final;

  /// Get property
    virtual std::string getPropertyS(std::string const& propertyName ) const override final;

  /// Get property
    virtual int  getPropertyI(std::string const& propertyName ) const override final;

 private:
    std::string _set_name{""};
    int _member{0};
    LHAPDF::PDF* _pdf{nullptr};
};

};  // namespace xfitter

