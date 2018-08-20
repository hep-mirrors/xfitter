#pragma once
#include "BasePdfDecomposition.h"
/**
  @class GRV_PionPdfDecomposition 

  @brief A class for GRV_Pion pdf decomposition
HACKS for s fraction

  Used for pi-
  Assumes that at starting scale:
    ubar=d
    dbar=u
    s=sbar
		u_sea=u   =f_u*sea
		d_sea=dbar=f_d*sea
		s_sea=s   =f_s*sea
		f_u=f_d
		f_u+f_d+f_s=1
		=> s=2*f_s/(1-f_s)*u
  Parametrised distributions are:
    v   := dval-uval=2*(d-u)
    qbar:=(u+dbar)/2=u
    g
  Therefore, transformations to physical basis:
    u=dbar=qbar
    d=ubar=qbar+v/2
    g=g
    s=sbar=qbar*2f_s/(1-f_s)
    others=0
  And sum rules for pi- are:
    \int_0^1 v dx=2
    \int_0^1 x*(4*qbar/(1-f_s)+v+g) dx=1
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
