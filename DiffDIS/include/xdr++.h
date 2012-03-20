/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef _XDR2CC_H_
#define _XDR2CC_H_

// #include <cstdio>
#include <cstdlib>
// #include <string.h>

extern "C" {
  #include <rpc/types.h>
  #include <rpc/xdr.h>
}

#define XDR_PTR ((XDR*)xdr)
/** 
  \brief C++ wrapper for basic XDR io.
  \internal
  \todo Add RWint16, RWint32.
  \endinternal
*/
// oooooooooooooooooooooooooooooooooooooooooooooo
class XDRio_t {
  void* xdr;
public:
  enum opmode {W,R} dir;
  
  XDRio_t(FILE* f, opmode mode) {
    xdr = NULL;
    xdr = malloc(sizeof(XDR));
    if(!xdr || !f) throw "XDRio_t creation failed";
    dir = mode;
    xdrstdio_create(XDR_PTR, f, (xdr_op)mode);
  }
  
  ~XDRio_t() {if(xdr) free(xdr); xdr=NULL;}
  
  bool isReading() {return dir;}
  
  //===================================
  bool RWchar(char* p) {return xdr_char(XDR_PTR,p);}
  bool RWuchar(unsigned char* p) {return xdr_u_char(XDR_PTR,p);}
  bool RWint(int* p) {return xdr_int(XDR_PTR,p);}
  bool RWdouble(double* p) {return xdr_double(XDR_PTR,p);}
  bool RWcstr(caddr_t *p, u_int maxlen) {return xdr_string(XDR_PTR,p,maxlen);}
  
  // bool RWchar(char* p);
  // bool RWuchar(unsigned char* p);
  // bool RWint(int* p);
  // bool RWdouble(double* p);
  // bool RWcstr(char**sp, unsigned int maxlen);
  // bool RWintv(int* p, unsigned int nelem);
  // bool RWdoublev(double* p, unsigned int nelem);
  
  //===================================
  bool RWintv(int* p, unsigned int nelem) {
    int k;
    for(k=0; k < nelem; k++,p++) if(!xdr_int(XDR_PTR,p)) return 0;
    return 1;
  }

  //===================================
  bool RWdoublev(double* p, unsigned int nelem) {
    int k;
    for(k=0; k < nelem; k++,p++) if(!xdr_double(XDR_PTR,p)) return 0;
    return 1;
  }
  
};
#undef XDR_PTR

#endif
