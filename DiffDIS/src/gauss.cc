/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#include <cmath>
#include <cstdio>
using namespace std;

#include "gauss.h"

const real_type Gauss_cnst=1.0e-15;
int gauss_error;

// ==========================
real_type Gauss16(
 RealFcn f,
 real_type a,
 real_type b,
 real_type eps)
{
#define GAPT 4
#define GAPT2 2*GAPT
  static real_type
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

  real_type c1,c2,u,s8,s16;
  real_type gsl,delta,aa,bb,y;
  //register
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
  if(fabs(s16-s8) <= eps*(1.0+fabs(s16))){
  gsl+=s16;
  aa=bb;
  goto l5;
 }
  y*=0.5;
  if(fabs(y)  >  delta) goto l2;

  printf("Gauss16 ... too high accuracy required [%g,%g]\n",a,b);
  gauss_error = 1;
  return gsl;
#undef GAPT
#undef GAPT2
}

// ===================================
real_type Gauss32(
 RealFcn f,
 real_type a,
 real_type b,
 real_type eps)
{
#define GAPT 8
#define GAPT2 2*GAPT
  static real_type
  w1[GAPT]={.189450610455068,
          .182603415044924, .169156519395003, .149595988816577,
          .124628971255534, .0951585116824928, .0622535239386479,
          .0271524594117541
  },
  x1[GAPT]={.0950125098376374,
          .281603550779259, .458016777657227, .617876244402644,
          .755404408355003, .865631202387832, .944575023073233,
          .989400934991650
  },

  w2[GAPT2]={.0965400885147278, .0956387200792749,
          .0938443990808046, .0911738786957639, .0876520930044038,
          .0833119242269468, .0781938957870703, .0723457941088485,
          .0658222227763618, .0586840934785355, .0509980592623762,
          .0428358980222267, .0342738629130214, .0253920653092621,
          .0162743947309057, .00701861000947010
  },
  x2[GAPT2]={.0483076656877383, .144471961582796,
          .239287362252137, .331868602282128, .421351276130635,
          .506899908932229, .587715757240762, .663044266930215,
          .732182118740290, .794483795967942, .849367613732570,
          .896321155766052, .934906075937740, .964762255587506,
          .985611511545268, .997263861849482
  };

  real_type gsl,delta,aa,bb,y,c1,c2,u,s8,s16;
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
  if(fabs(s16-s8) <= eps*(1.0+fabs(s16))){
  gsl+=s16;
  aa=bb;
  goto l5;
 }
  y*=0.5;
  if(fabs(y)  >  delta) goto l2;

  printf("Gauss32 ... too high accuracy required\n");
  gauss_error = 1;
  return gsl;
#undef GAPT
#undef GAPT2
}
