
#pragma once

#include "BasePdfDecomposition.h"
#include "BasePdfParam.h"
#include <memory>

namespace xfitter
{
  /**
     @class UvDvUbarDbarS_PdfDecomposition 

     @brief A class for UvDvUbarDbarS_ pdf decomposition

     @version 0.1
     @date 2018-07-11
  */

  class UvDvUbarDbarS : public BasePdfDecomposition
  {
  public:
    /// Default constructor. Name is the PDF name
    UvDvUbarDbarS();

    /// Optional initialization at the first call
    virtual void initAtStart(const std::string& pars) override final;

    /// Compute sum-rules
    virtual void initAtIteration() override final;
    
    /// Compute PDF in a physical base in LHAPDF format at the initial scale
    virtual std::function<std::map<int,double>(const double& x)> f0() const override final;


  private:
    /// Get parameter values from the minimizer
    std::unique_ptr<double[]> getParValues(BasePdfParam const* param) const;

    /// Get valence
    double valence(double x, std::string const& name, double sum) const;

    /// Get sea
    double sea(double x, std::string const& name) const;
    
    /// Get uv, apply sum-rule
    double uv(double x) const;

    /// Get dv, apply sum-rule
    double dv(double x) const;    

    /// Get dbar
    double dbar(double x) const;

    /// Get ubar
    double ubar(double x) const;

    /// Get s
    double s(double x) const;

    /// Get g
    double g(double x) const;

  private:
    /// sum-rule fixed parameters
    double _uSum,  _dSum, _gSum;
    
  };
}
