#include "BaseMinimizer.h"

namespace xfitter {

  void BaseMinimizer::addParameterBlock(int Npar, double const* pars, std::string const* names,  unsigned int const* flags , double const* const* bounds  ) {
    for (int i =0 ; i<Npar; i++) {
      addParameter(pars[i], names[i],flags[i],bounds[i]);
    }
  }

  void BaseMinimizer::addParameter(double par, std::string const &name,  unsigned int flag, double  const* bounds  ) {
    _allParameters.push_back(par);
    _allParameterNames.push_back(name);
    if ( flag != isFixedPar ) {
      _fitParameters.push_back( &_allParameters.back() );
    }
  }
}

