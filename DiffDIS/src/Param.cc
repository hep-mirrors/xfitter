/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#include "Param.h"

// string ParamVec_t::BadChars(Xstring::white_space+"=&"); //--- this crashes on RedHat
string ParamVec_t::BadChars(" \t\r=&");
char Param_t::QuoteChar[2];
