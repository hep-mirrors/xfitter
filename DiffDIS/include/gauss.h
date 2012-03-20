/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef _WS_GAUSS_H_
#define _WS_GAUSS_H_

/** @addtogroup gr_mathut
 *  @{
 */

#include "types.h"
extern int gauss_error;

/**
\brief 16-point Gauss quadrature.
*/
real_type Gauss16(RealFcn f, real_type a, real_type b, real_type eps=1e-6);
/**
\brief 32-point Gauss quadrature.
*/
real_type Gauss32(RealFcn f, real_type a, real_type b, real_type eps=1e-6);

///@}

#endif
