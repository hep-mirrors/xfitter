/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef QCDPARAMS_H_
#define QCDPARAMS_H_

#include "types.h"

#define N_SFTYPES 6
enum SFTYPE {FF,FG,GF,GG,PLUS,MINUS};
/// @cond NIEWAZNE
extern const real_type eq2[], c_F, c_A;
extern const char* SFTname[];
/// @endcond

#endif
