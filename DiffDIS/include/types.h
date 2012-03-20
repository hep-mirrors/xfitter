/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef TYPES_H_
#define TYPES_H_

/** @addtogroup gr_mathut
 *  @{
 */

// typedef long double real_type;
typedef double real_type;

typedef real_type*  RealPtr_t;
typedef RealPtr_t* RealPtr2_t; //--- 2 dim array
typedef const real_type*  ConstRealPtr_t;
typedef const RealPtr_t* ConstRealPtr2_t;
typedef real_type (*RealFcn)(real_type);

///@}

#endif
