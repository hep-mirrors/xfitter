/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

/* Lambert W function, -1 branch. 

   K M Briggs 2001 Apr 09
   Keith dot Briggs at bt dot com

   Returned value W(z) satisfies W(z)*exp(W(z))=z
   Valid input range: -1/e <= z    < 0
   Output range:      -1   >= W(z) > -infinity

   See LambertW.c for principal branch.
   Together these form the only purely real branches.

   Test: 
     gcc -DTESTW LambertW1.c -o LambertW1 -lm && LambertW1
   Library:
     gcc -O3 -c LambertW1.c 
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>

// double LambertW1(const double z) {
double LambertW1(double z) {
  int i; 
  const double eps=4.0e-16, em1=0.3678794411714423215955237701614608;
  // >> --- ws
  // double p,e,t,w,l1,l2;
  double p=1;
  double e,t,w,l1,l2;
  // <<
  if (z<-em1 || z>=0.0 /*|| isinf(z) || isnan(z)*/) {
    // fprintf(stderr,"LambertW1: bad argument %g, exiting.\n",z); 
    // exit(1);
    char ts[64];
    sprintf(ts,"LambertW1: bad argument %g",z);
    throw ts;
  }
  /* initial approx for iteration... */
  if (z<-1e-6) { /* series about -1/e */
    p=-sqrt(2.0*(2.7182818284590452353602874713526625*z+1.0));
    w=-1.0+p*(1.0+p*(-0.333333333333333333333+p*0.152777777777777777777777)); 
  } else { /* asymptotic near zero */
    l1=log(-z);
    l2=log(-l1);
    w=l1-l2+l2/l1;
  }
  if (fabs(p)<1e-4) return w;
  for (i=0; i<10; i++) { /* Halley iteration */
    e=exp(w); 
    t=w*e-z;
    p=w+1.0;
    t/=e*p-0.5*(p+1.0)*t/p; 
    w-=t;
    if (fabs(t)<eps*(1.0+fabs(w))) return w; /* rel-abs error */
  }
  /* should never get here */
  // fprintf(stderr,"LambertW1: No convergence at z=%g, exiting.\n",z); 
  // exit(1);
  char ts[64];
  sprintf(ts,"LambertW1: No convergence at z = %g",z);
  throw ts;
}

#ifdef TEST_LAMBERTW1
/* test program...  */
int main() {
  int i;
  double z,w,err;
  for (i=0; i<100; i++) {
    z=i/100.0-0.3678794411714423215955; 
    if (z>=0.0) break;
    w=LambertW1(z); 
    err=exp(w)-z/w;
    printf("W(%8.4f)=%22.16f, check: exp(W(z))-z/W(z)=%e\n",z,w,err);
  }
  z=-1.0e-2; 
  for (i=0; i<10; i++) {
    z/=10.0; 
    w=LambertW1(z); 
    err=exp(w)-z/w;
    printf("W(%8.4g)=%22.16f, check: exp(W(z))-z/W(z)=%e\n",z,w,err);
  }
  return 0;
}
#endif
