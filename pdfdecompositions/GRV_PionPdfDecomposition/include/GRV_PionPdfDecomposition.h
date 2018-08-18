#pragma once
#include "BasePdfDecomposition.h"
/**
  @class GRV_PionPdfDecomposition 

  @brief A class for GRV_Pion pdf decomposition

  Used for pi-
  Assumes that at starting scale:
    ubar=d
    dbar=u
    s=sbar=0
  Parametrised distributions are:
    v   := dval-uval=2*(d-u)
    qbar:=(u+dbar)/2=u
    g
  Therefore, transformations to physical basis:
    u=dbar=qbar
    d=ubar=qbar+v/2
    g=g
    others=0
  And sum rules for pi- are:
    \int_0^1 v dx=2
    \int_0^1 x*(4*qbar+v+g) dx=1
  @version 0.2
  @date 2018-08-14
  */
namespace xfitter{
class GRV_PionPdfDecomposition:public BasePdfDecomposition{
	public:
		GRV_PionPdfDecomposition();
		GRV_PionPdfDecomposition(const std::string& inName);
		~GRV_PionPdfDecomposition();
		virtual void initAtStart(const std::string & pars) override final;
		virtual void initAtIteration()override final;
		virtual std::function<std::map<int,double>(const double& x)>f0()const override final; 
	private:
		BasePdfParam*par_v,*par_qbar,*par_g;
	};
}
