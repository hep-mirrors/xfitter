/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

/* \file dilog.h
  \brief Spence, dilogarithm, S2
*/
#include "dilog.h"
const real_type pi2o6 = M_PI*M_PI/6;

/**
  \brief Re Spence(x)
  i.e. for x < -1 returns Real part only
  
  \f[
  {\rm Spence}(x) = \int\limits_0^x \! dz\; \frac{\ln(1+z)}z
  \f]
*/
//=====================================================
real_type Spence(real_type x){
  const real_type eps=1e-16;
  // printf("Spence for: %lf = ", x);
  
  // --- special arg. values: -1, -0.5, 1
  if(fabs(x+1) < eps) return -pi2o6;
  if(fabs(x+0.5) < eps) return 0.5*(M_LN2*M_LN2 - pi2o6);
  if(fabs(x-1) < eps) return 0.5*pi2o6;

  real_type sum=0, w, xneg=-x;
  int n=1;
  if(x < -0.5) goto D1;
  if(x < -0.3) goto D2;
  if(x < 0.0)  goto NN;
  if(x > 1.0)  goto RR;
  if(x > 0.7)  goto NN;

  // --- Abramovitz, Stegun (27.7.2)
  do {
    sum += w = pow(xneg,n)/(n*n);
    n++;
  } while(fabs(w) >= eps);
  return -sum;

  D1:  return log(-x)*log(fabs(1+x)) -pi2o6 - Spence(-x-1);
  D2:  return 0.5*Spence(-x*x) - Spence(-x);
  NN:  return 0.5*pow(log(1+x),2) - Spence(-x/(1+x));
  RR:  return 0.5*pow(log(x),2) + pi2o6 - Spence(1.0/x);
}

/// dilog(x) = -Spence(x-1)
/**
  \brief Re dilog(x)
  i.e. for x < 0 returns Real part only
  
  \f[
  {\rm dilog}(x) = -\int\limits_1^x \! dz\; \frac{\ln z}{z-1}
  \f]
*/
//=====================================================
real_type dilog(real_type x){
  return -Spence(x-1);
}

//=====================================================
real_type S2(real_type x){
  real_type lnx=log(x);
  return 2.0*Spence(x)-pi2o6-lnx*(2.0*log(1.0+x) - 0.5*lnx);
}
