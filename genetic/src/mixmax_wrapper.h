// definitions for the MIXMAX generator

#ifndef __LP64__
typedef uint64_t myuint;
#else
typedef unsigned long long int myuint;
#endif

extern "C" {
  myuint* rng_alloc(void);
  void seedvector(myuint*, int);
  long double get_next_float(myuint*);
  void rng_free(myuint*);
  myuint* skip_1M(myuint*);
}

// functions
void ini_rnd(void);
int get_uniform_int(int min, int max);
double get_uniform_double(double min, double max);
double get_gaus_double(double mean, double sigma);
void finalize_rnd(void);

