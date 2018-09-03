#pragma once

#include  "BaseEvolution.h"


#include <vector>
#include <memory>
#include <functional>


/**
  @class' EvolutionQCDNUM

  @brief A wrapper class for QCDNUM evolution 

  @version 0.2
  @date 2018-09-29
  */

namespace xfitter
{

  class EvolutionQCDNUM: public BaseEvolution
  {
  public:
    EvolutionQCDNUM(const char*name):BaseEvolution{name}{};

    virtual void initAtStart()override final;
    virtual void initAtIteration()override final;
    virtual void initAtParameterChange()override final;
    virtual std::function<std::map<int,double>(double const& x, double const& Q)> xfxQMap()  override final;
    virtual std::function<void(double const& x, double const& Q, double* pdfs)> xfxQArray() override final;
    virtual std::function<double(int const& i, double const& x, double const& Q)> xfxQDouble() override final;
    virtual std::function<double(double const& Q)> AlphaQCD() override final;
//  protected:
//    virtual int parseOptions(){ return 0;};
  private:
    /// PDFs called outside boundaries check:
    int _icheck{0};

    /// Number of external PDFs
    int _nExt{0};

    /// PDF type (1 -- unpolorized internal)
    int _itype{1};

    /// Spline order for interpolation
    int _splineOrder{2};

    /// Read or not QCDNUM pre-stored tables (makes init faster)
    int _readTables{0};
    
  };
};


