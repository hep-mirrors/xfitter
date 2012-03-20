/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

/* \file mathut.h
  \brief math utilities
  
  Constants, types, classes, functions, ...
*/
#ifndef WS_MATH_UT
#define WS_MATH_UT

/** @defgroup gr_mathut Mathematical utilities
 *  Constants, types, functions and classes
 *  @{
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;

#ifndef M_PI
  #define M_PI        3.1415926535897932384626433832795
#endif
#ifndef M_PI_SQ
  #define M_PI_SQ     9.86960440108935861883449099987615
#endif

#include "types.h"

//***********************************************************************
template <typename T> int sign(T x)
{return (x < 0) ? -1 : 1;}


//***********************************************************************
template <class T> int outrange(const T& x, const T& xmin, const T& xmax)
{return (x < xmin) ? -1 : ((x >= xmax)? 1 : 0);}

/// Interval of \c real_type
//-----------------------------------------------------------------------
class TLimits {
// public:
  real_type lo, hi;
protected:
  int started;

public:
  TLimits(){ started = 0;}
  void Reset(){started = 0;}
  void Reset(real_type xl, real_type xh) {
    lo = xl;
    hi = xh;
    started = 1;
  }

  TLimits(real_type xl, real_type xh) {
    Reset(xl, xh);
  }

  real_type Lo() const {
    if(!started) {
      cerr << "TLimits class ERROR - not started in Lo.\nPress <Enter> to close." << endl;
      cin.get();
      exit(101);
    }
    return lo;
  }

  real_type Hi() const {
    if(!started) {
      cerr << "TLimits class ERROR - not started in Hi.\nPress <Enter> to close." << endl;
      cin.get();
      exit(101);
    }
    return hi;
  }

  void IncrLo(real_type x){
    if(x > lo) lo = x;
  }
  void DecrHi(real_type x){
    if(x < hi) hi = x;
  }

  void Intersect(real_type xl, real_type xh){
    if(!started) return;
    IncrLo(xl);
    DecrHi(xh);
  }

  void Intersect(const TLimits& lims){
    if(!started) return;
    IncrLo(lims.lo);
    DecrHi(lims.hi);
  }

  TLimits Scale(real_type s){
    if(!started) return *this;
    if(s >= 0.0) {lo = s*lo; hi = s*hi;}
    else {real_type lo1 = s*hi; hi = s*lo; lo = lo1;}
    return *this;
  }

  void Expand(real_type x){
    if(started){
      if(x < lo) lo = x;
      else if(x > hi) hi = x;
    }
    else {lo = hi = x; started = 1;}
  }
/*
  void Expand(real_type xl, real_type xh){
    Expand(xl);
    Expand(xh);
  }
*/
  bool InRange(real_type x) const {
    return started && (lo < x) && (x < hi);
  }

  bool InRangeR(real_type x) const {
    return started && (lo < x) && (x <= hi);
  }

  bool InRangeL(real_type x) const {
    return started && (lo <= x) && (x < hi);
  }

  bool InRangeLR(real_type x) const {
    return started && (lo <= x) && (x <= hi);
  }

  friend ostream& operator << (ostream& ops, const TLimits& rr) {
    return ops << "[" << rr.Lo() << ", "
        << rr.Hi() << "] ";//<< setw(DBLWIDTH)
  }

};

/// A sequence of linearly or logarithmicaly spaced values
// ooooooooooooooooooooooooooooooo
class TTableVar {
  real_type step, start, end, v;
  int logstep;

public:
  TTableVar(real_type s, real_type e, int np, int logst=0) {
    start = s;
    end = e;
    logstep = logst? 1 : 0;
    if(np > 1) step = logstep? log(end/start)/(np-1) : (end-start)/(np-1);
    else step = logstep;
  }

  TTableVar(const TLimits& lims, int np, int logst=0) {
    *this = TTableVar(lims.Lo(), lims.Hi(), np, logst);
  }

  /*
  real_type Step(int np){
    return step = (end-start)/(np-1);
  }
  */
  real_type Val(int ind) {
    return v = logstep? start*exp(step*ind) : (start + ind*step);
  }
};


//********************************************************
inline real_type cosbar(real_type theta) {
  real_type sh = sin(theta*0.5);
  return 2.0*sh*sh;
}

#include "gauss.h"
#include "dilog.h"

//--- in LambertW1.cc
double LambertW1(double z);

//--- in tblint.cc
double TblInt1(double xh, double yi[], int np);
double TblInt2(double xi[], double yi[], int np);

#include "tblipol.h"
#include "MDipol.h"

/** @} */ // end of gr_mathut

#endif
