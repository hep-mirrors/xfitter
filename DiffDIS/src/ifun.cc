/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#include "ifun.h"

using namespace std;

//================================
double ifun_t::Gauss16(double a, double b, double acc) {
static const double Gauss_cnst=1.0e-15;
#define GAPT 4
#define GAPT2 2*GAPT
  static double
  w1[GAPT]={
    .362683783378362, .313706645877887,
    .222381034453374, .101228536290376
  },
  x1[GAPT]={
    .183434642495650, .525532409916329,
    .796666477413627, .960289856497536
  },
  w2[GAPT2]={.189450610455068,
          .182603415044924, .169156519395003, .149595988816577,
          .124628971255534, .0951585116824928, .0622535239386479,
          .0271524594117541
  },
  x2[GAPT2]={.0950125098376374,
          .281603550779259, .458016777657227, .617876244402644,
          .755404408355003, .865631202387832, .944575023073233,
          .989400934991650
  };

  double c1,c2,u,s8,s16;
  double gsl,delta,aa,bb,y;
  int i;

  gauss_error = 0;
  delta=Gauss_cnst*fabs(a-b);
  gsl=0.0;
  aa=a;
l5: y=b-aa;
  if(fabs(y)  <=  delta) return gsl;
l2: bb=aa+y;
  c1=0.5*(aa+bb);
  c2=c1-aa;
  s8=0.0;
  s16=0.0;
  for(i=0;i<GAPT;i++) {
	  u=x1[i]*c2;
	  s8+= w1[i]*(f(c1+u)+f(c1-u));
  }
  for(i=0;i<GAPT2;i++){
	  u=x2[i]*c2;
	  s16+= w2[i]*(f(c1+u)+f(c1-u));
  }
  s8*=c2;
  s16*=c2;
  if(fabs(s16-s8) <= acc*(1.0+fabs(s16))){
  gsl+=s16;
  aa=bb;
  goto l5;
 }
  y*=0.5;
  if(fabs(y)  >  delta) goto l2;

  fprintf(stderr, "Gauss16 ... too high accuracy required [%g,%g]\n",a,b);
  gauss_error = 1;
  return gsl;
#undef GAPT
#undef GAPT2
}
