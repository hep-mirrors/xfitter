/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef _BASIC_DEFS_H_
#define _BASIC_DEFS_H_

#include "genut.h"
#include "dbgtools.h"

#define MAX_FLAVORS 6
// --- -6,...,7 with a space for q_sum
#define N_FLAV_INDICES (2*MAX_FLAVORS+2)
enum Parton_t {gluon, d_quark, u_quark, s_quark, c_quark, b_quark, t_quark, q_sum};
// enum Parton_t {gluon, d_quark, u_quark, s_quark, c_quark, b_quark, t_quark};

#define MAX_QCD_ORDER 1

// --- default values
#define MASS_C 1.5
#define MASS_B 4.5
#define MASS_T 174.0
#define MASS_Z0 91.1876
#define ALPHAS_Z0 0.118

#endif
