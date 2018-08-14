#pragma once

#include  "BaseEvolution.h"


#include <vector>
#include <memory>
#include <functional>


/**
  @class' EvolutionQCDNUM

  @brief A wrapper class for QCDNUM evolution 

  @version 0.1
  @date 2018-08-14
  */

namespace xfitter
{

  class EvolutionQCDNUM: public BaseEvolution
  {
  public:
    /// Empty constructor (needed for the dynamic loading)
    EvolutionQCDNUM():  BaseEvolution{"QCDNUM",nullptr}{};
    
  public:
    virtual std::string getEvolutionName() const { return  "QCDNUM" ;};
    virtual void initAtStart() override final;
    virtual void initAtIteration() override final;
    virtual std::function<std::map<int,double>(double const& x, double const& Q)> xfxQMap()  override final;
    virtual std::function<void(double const& x, double const& Q, double* pdfs)> xfxQArray() override final;
    virtual std::function<double(int const& i, double const& x, double const& Q)> xfxQDouble() override final;
    virtual std::function<double(double const& Q)> AlphaQCD() override final;
  protected:
    virtual int parseOptions(){ return 0;};
  };
};


