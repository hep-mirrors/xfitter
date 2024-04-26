#pragma once
#include "BasePdfDecomposition.h"

#include "spline.h"


/**
  @class Q0decomposition

  @brief A class for pdf decomposition at a starting scale from a text file

  @version 0.1
  @date 2024-04-25
**/
namespace xfitter
{
  class Q0decomposition : public BasePdfDecomposition
  {
  public:
    Q0decomposition(const char*name);
    ~Q0decomposition();
    virtual const char*getClassName()const override final;
    virtual void atStart()override final;
    virtual std::map<int,double>xfxMap(double x)const override final;
  private:
    std::map <int, tk::spline> _PDFsAtQ0;
  };
}
