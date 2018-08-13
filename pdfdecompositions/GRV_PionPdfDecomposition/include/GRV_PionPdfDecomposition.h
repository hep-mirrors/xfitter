
#pragma once

#include "BasePdfDecomposition.h"

/**
  @class GRV_PionPdfDecomposition 

  @brief A class for GRV_Pion pdf decomposition

  @version 0.1
  @date 2018-08-07
  */

/*
Used for pi-
Assumes that at starting scale:
	ubar=d
	dbar=u
	s=sbar=0
Parametrised distributions are:
	v   := dval-uval
	qbar:=(ubar+dbar)/2
	g
Therefore, transformations to physical basis:
	u=dbar=qbar-v/4
	d=ubar=qbar+v/4
	g=g
	others=0
And sum rules for pi- are:
	\int_0^1 v dx=2
	\int_0^1 x*(4*qbar+g) dx=1
*/

namespace xfitter {
class ParameterisationWrapper;
class GRV_PionPdfDecomposition:public BasePdfDecomposition{
	public:
		GRV_PionPdfDecomposition();
		 /// Default constructor. Name is the PDF name
		GRV_PionPdfDecomposition(const std::string& inName);
		~GRV_PionPdfDecomposition();
		/// Optional initialization at the first call
		virtual void initAtStart(const std::string & pars) override final;
		/// Compute PDF in a physical base in LHAPDF format for given x and Q
		virtual std::function<std::map<int,double>(const double& x)>f0()const override final; 
	private:
		ParameterisationWrapper*par_v,*par_qbar,*par_g;
	};
}
