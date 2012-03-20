/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef NO_XDR
#include "GrdHeader.h"

//==========================================================
bool GrdHeader_t::RW(XDRio_t& xdr) {
  bool OK=1;
  //fread(&version,
  //return 1;
  //xdr_array(&xdrs, &Vp, &len, NVALS, sizeof(V[0]), (xdrproc_t)xdr_double);
  char *sp = lbl, *sp2=info;
  OK = OK
    && xdr.RWuchar(&Ver)
    && xdr.RWuchar(&SubVer)
    && xdr.RWcstr(&sp, 63)
    && xdr.RWcstr(&sp2, 127)
    && xdr.RWint(&QCDClib_ver)
    && xdr.RWint(&Ndistr)
    && xdr.RWint(&x_exp)
    && xdr.RWint(&PiOrd)
  ;

  return OK;
}

#endif
