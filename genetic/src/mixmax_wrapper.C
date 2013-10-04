#include "mixmax_wrapper.h"
#include <math.h>

// keeps a current position
myuint* V;

void ini_rnd(void)
{
  // initialize MIXMAX generator
  V = rng_alloc();
  seedvector(V, 1);

  // TODO: use the std function?
  for (int i = 0; i < 1000; ++i) {
    (void)get_next_float(V);
  }
  //skip_1M(V); // use this instead

}

int get_uniform_int(int min, int max)
{
  // the return value is in [min, max]

  // shouldn't it be several idle calls here?
  return ((int)(get_next_float(V)*(max - min + 1) + min));
}

double get_uniform_double(double min, double max)
{
  // the return value is in [min, max]

  // shouldn't it be several idle calls here?
  return (((double)get_next_float(V))*(max - min) + min);
}

double get_gaus_double(double mean, double sigma)
{
  // Use Box-Muller algorithm
  // caution: be aware of improper tails

  static bool even = true;

  static double z1 = 0.;
  static double z2 = 0.;

  if (even) {
    
    double x = get_uniform_double(-1.,1.);
    double y = get_uniform_double(-1.,1.);
    double s = x*x + y*y;
    while (s > 1. || s == 0.) {
      x = get_uniform_double(-1.,1.);
      y = get_uniform_double(-1.,1.);
      s = x*x + y*y;
    }

    double logs = sqrt(-2.*log(s) / s);
    z1 = x*logs;
    z2 = y*logs;

    even = false;

    return (mean + sigma*z1);
  } else {

    even = true;
    return (mean + sigma*z2);
  }

}

void finalize_rnd(void)
{
  rng_free(V);
}
