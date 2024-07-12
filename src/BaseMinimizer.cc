#include "BaseMinimizer.h"
#include "xfitter_pars.h"
#include <memory>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace xfitter {

  void BaseMinimizer::addParameterBlock(int Npar, double const* pars
				   , std::string const* names
				   , double const* steps   
				   , double const* const* bounds 
				   , double const* const* priors  )
  {
    for (int i =0 ; i<Npar; i++) {
      addParameter(pars[i], names[i],steps[i],bounds[i],priors[i]);
    }
  }

  void BaseMinimizer::addParameter(double par, std::string const &name, double step, double const* bounds , double  const* priors  )
  {

    // names for minimized parameters only
    if ( step > 0) {
      _allParameterNames.push_back(name);
      _bounds.push_back(bounds);
      _priors.push_back(priors);    
    }

    // store it on the global map too. Will replace pointer if already present. 
    //    std::unique_ptr<double[]> parval( new double );

    if ( XFITTER_PARS::gParameters.find(name) !=  XFITTER_PARS::gParameters.end() ) 
      *XFITTER_PARS::gParameters[name] = par;
    else {    
      double*parval = new double;
      *parval = par;    
      XFITTER_PARS::gParameters[name] = parval;
    }
  }

  void BaseMinimizer::CopyStateFromMinimizer(const BaseMinimizer* origin){
    for (int ipar=0; ipar<origin->_allParameterNames.size(); ipar += 1) {
      std::string name = origin->_allParameterNames[ipar];
      double value = *XFITTER_PARS::gParameters[name];
      const double* bounds = origin->_bounds[ipar];
      const double* priors = origin->_priors[ipar];
      addParameter(value, name, fabs(value)*1.e-2, bounds, priors);
    }
  }

  double** BaseMinimizer::getPars() const
  {
    double**   out =  new double*[getNpars()] ;
    for (size_t i = 0; i<_allParameterNames.size(); i++) {      
      out[i] = XFITTER_PARS::gParameters.at( _allParameterNames[i] );
    }
    return out;
  }

  void BaseMinimizer::setPars(double const* pars) const {
    std::cout << std::setw(5) << "NO" << std::setw(30) << "NAME" << std::setw(15) << "VALUE" << std::endl;
    for (size_t i = 0; i<_allParameterNames.size(); i++) {
      std::cout << std::setw(5) << i << std::setw(30)  << _allParameterNames[i] << std::setw(15) <<  pars[i] << std::endl;
      //std::cout << i << " parval = " << pars[i] << std::endl;
      *XFITTER_PARS::gParameters.at( _allParameterNames[i]) = pars[i];
    }
  }
}

