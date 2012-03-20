/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef TBLIPOL_H_
#define TBLIPOL_H_

/** @addtogroup gr_mathut
 *  @{
 */

#include "types.h"
int FindIndex(const real_type xx[], int xx_size, real_type x);
int FindIndexR(const real_type xx[], int xx_size, real_type x);
real_type Ipol1(real_type x, const real_type xx[], const real_type yy[], int npt, int pol_deg=2);
// typedef const real_type* crp_t;
real_type Ipol2(real_type x, real_type y,
  // const real_type xx[], const real_type yy[],  const crp_t* zz, int nptx, int npty
  const real_type xx[], const real_type yy[],  ConstRealPtr2_t zz, int nptx, int npty
  , int pol_degx=2, int pol_degy=2);

///@}

#endif
