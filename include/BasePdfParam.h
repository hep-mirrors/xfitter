#pragma once


#include <string>
#include <yaml-cpp/yaml.h>

/**
  @class BasePdfParam

  @brief Base class for PDF parameterisations

  Represents a function xf(x) that depends on x and some external parameters --- xf(x,A,B,C,...)
  Has methods to
  *Evaluate the function at given x and current (globally stored) parameters
  *Calculate n-th moment, e.g. \int_0^1 x^n*xf(x) dx
  *Rescale some parameters so that the n-th moment has a given value

  Parameterisations keeps pointers to its parameters.

  Parameterisations are initialised by parsing instance-specific YAML nodes
  Concrete implementations of BasePdfParam must be able to get their parameters and any/all additional options

  @version 0.2
  @date 2018-08-13
  */

namespace xfitter{
class BasePdfParam{
public:
  BasePdfParam(const std::string&instance_name):_name(instance_name){}
  virtual~BasePdfParam();
  void               setNPar(unsigned int N){Npars=N;}
  const unsigned int getNPar()const{return Npars;}
  //!Evaluate xf(x) at given x with current parameters, pure virtual
  virtual double operator()(double x)const=0;
  //!Calculate n-th moment, e.g. \int_0^1 x^n*xf(x) dx
  //!Note that we parameterise xf(x), not f(x)
  //!Therefore, moment(-1) should be used for valence sums, and moment(0) for momentum sum
  virtual double  moment(int nMoment=-1)const;
  //!Rescale some parameters so that the n-th moment has a given value
  //!A typical implementation will probably do this by setting *pars[0]
  virtual void setMoment(int nMoment,double value);
  //!Get name of the instance
  const std::string getName()const{return _name;}
  virtual void atStart();
  virtual void atIteration();
protected:
  //!Unique name of instance
  const std::string _name;
  //!Array of pointers to some global locations where minimization parameters are stored
  double**pars{nullptr};
  //!Number of parameters, which is also the size of the array **pars defined above
  unsigned int Npars{0};
};
}
