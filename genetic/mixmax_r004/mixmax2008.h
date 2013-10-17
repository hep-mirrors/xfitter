#include <stdint.h>

#define BITS  61

/* Mersenne Numbers */
#define M3    7
#define M5    31
#define M13   8191
#define M31   2147483647
#define M61   2305843009213693951
#define M(B)  M##B
#define xM(B) M(B)
#define BASE  xM(BITS)

#ifndef N
#define N 60 //128
/* For N=32, will generate (and print) a billion numbers in under 4min, on Intel Xeon 2.8G . 
   Printing, I think is the more expensive part! 
   Since the algorithm is linear in N, the cost per number is roughly independent of N.
   */
#endif

//#ifndef __LP64__
#ifdef __ILP32__
typedef uint64_t myuint;
#else
typedef unsigned long long int myuint;
#endif

myuint *rng_alloc();
int rng_free(myuint *X);

int seedvector(myuint *X, unsigned int seed);
int seed_lcg(myuint *X);

inline void iterate(myuint *X);
inline void iterate_unrolled(myuint *X);
myuint get_next(myuint *X);
long double get_next_float(myuint *X);

myuint *skip_1M(myuint *X);
myuint *skip_1G(myuint *X);


#define STR(a) X ## a
#define XSTR(x) STR(x)
#define NS N
#define ALL(a)	XSTR(NS)(a)
/*  
 *  These macros are there to allow automatic code generation for unrolled loops
 *  up to size N=256.  N should be power of 2 if using these macro. 
 */

#define IT(a) iterate(a)
/* iterate or iterate_unrolled (version using macros seems to be three times faster)
 * at N=8, macro version computes 10^8 iterations in 0m8.376s
 * while normal loopy version in 0m25.647s
 * however, optimizer with -funroll-loops (or -O3) seems to be fastest at 0m6.693s
 */


#define SUFF(i) i##ULL
#define ZERO SUFF(0)
#define ONE SUFF(1)
#define SEEDVAL SUFF(1)
//2^63=9223372036854775808, 2^62=4611686018427387904, 2^60=1152921504606846976, 2^58=288230376151711744


#define MODULAR 1

#ifdef MODULAR
//#define MERSBASE SUFF(2305843009213693951) // 2^61-1 = 2305843009213693951
#define MERSBASE BASE // 2^61-1 = 2305843009213693951
#define MOD_MERSENNE(k) ((k) - ((k+1) >> BITS)*MERSBASE ) 
#define INV_MERSBASE .000000000000000000433680868994201773791060216479542685926876L // that is 1/(2^61-1), for numbers in [0,1), if numbers in [0,1] are desired, use 1/(2^61-2), but remember, this only matter only once in a femtillion.
#else
#define MOD_MERSENNE(k) (k)
#define INV_MERSBASE .000000000000000000054210108624275221700372640043497085571289L // that is 1/(2^64)
#endif




#define FLOAT_OUTPUT 1
#ifndef FLOAT_OUTPUT
#define OUT(x) printf("%20llu,\t", x);
#else
#define OUT(x) printf("%1.20LF\t", (long double)x* INV_MERSBASE); // 20 digits of which 18 or 19 are good
#endif
