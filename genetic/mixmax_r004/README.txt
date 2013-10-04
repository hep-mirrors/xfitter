README

Welcome to the new, computationally super-efficient O(N) implementation of MIXMAX, the Yerevan random number generator!

MIXMAX is a matrix-recursive random number generator introduced in 1986 by my parents, George Savvidy and Natalia Ter-Arutyunyan-Savvidy.

[1]   On the Monte Carlo simulation of physical systems
      J.Comput.Phys. 97, 566 (1991), 
      http://dx.doi.org/10.1016/0021-9991(91)90015-D
      http://ccdb5fs.kek.jp/cgi-bin/img_index?8607219

and, the implementation, with N. Akopov:

[2] Matrix generator of pseudorandom numbers
    J.Comput.Phys.97, 573 (1991)
    http://dx.doi.org/10.1016/0021-9991(91)90016-E
    http://ccdb5fs.kek.jp/cgi-bin/img/allpdf?198607220

WHY USE MIXMAX?

If you are doing Monte-Carlo in high-dimensional space, you should consider using MIXMAX. Examples are particle physics, lattice gauge theory, statistical physics or whenever the number of random numbers to simulate each event, trajectory, MCMC or Metropolis update and so on is large. How large? If it is any larger than the state vector of the generator you must not use it. To give an idea, Luscher's widely used generator has a state vector of size 24. The various Mersenne twisters are larger, but not recommended for physics due to very long correlation lengths in their sequences, as explained below in the THEORY section. Unlike generators designed by computer scientists which tend to emphasize equi-distribution and randomness in terms of bits, MIXMAX has strong theoretical guarantees in terms of its floating point output. The default dimension used in this implementation is N=60, but we invite the user to increase it as necessary if his/her Monte-Carlo dimension is larger.

WHAT IS MIXMAX?

The random number generator was an outgrowth of research on dynamical systems, namely Yang-Mills classical mechanics. It was noted in the first paper that if a dynamical system possessed the property of Kolmogorov's mixing, then it could be used to generate pseudo-random numbers. A specific system with discrete time was proposed [1], by means of a linear automorphism of a hypercube, x_{i+1} = A . x_i mod 1 where A is a matrix with integer entries. For mixing, it is required that the determinant is unit, and that all eigenvalues are different by absolute value from 1. A specific matrix  [2] with these properties is

2 3 4 5 É    N  1
1 2 3 4 É   N-1 1
  ...
1 1 1 1 ... 2  2  1       <- this is not a typo, the second entry from the right is 2, not 3!
1 1 1 1 ... 1  2  1
1 1 1 1 ... 1  1  1

The original generator was thus based on a recursion with real numbers on [0,1), generating directly in floating point. The present implementation allows to do this too, but by default applies the same  recursion to vectors with integer entries a_i modulo a Mersenne prime, here p=2^61-1. As a dynamical system it is close, in a precise sense, to the continuous one, we identify   x_i  <==> a_i / p and output x_i if floating point numbers are desired. 


INSTALLATION

On Linux, Unix, and Mac OSX systems, unzip the archive and change into the directory:

unzip mixmax_release_NNN.zip
cd mixmax_release_NNN

The implementation of the generator is in the file mixmax.c and there are various example driver programs included. Type make, and hopefully you will have an executable called mixmax2008 which can be run:

make
./mixmax2008

At this point, it will ask you to enter the number of floating point numbers to produce, and then print them one on a line, in 20 decimal digit precision which is a bit more than the actual usable precision of 18 digits in this implementation of the generator. In actual use, you are permitted to request floating point with get_next_float() and integer numbers with get_next() in any order you may need. If you need 32 bit or some other integers, it is recommended to bit shift the integer at hand by 29 bits to the right.


PORTABILITY

The implementation uses the 64 bit unsigned integer which should be available on most modern systems, including some 32 bit systems. The actual arithmetic is modular in the Mersenne prime base 2^61-1, and is close in a precise sense to the dynamical system originally described in those two papers. Having said that, the actual period of the cycle of this generator is expected to be of order p^N but the precise value is unknown for arbitrary N. For N=16 it has been recently verified that the period is maximal, p^N -1. With modular arithmetic it is guaranteed that there are no bad seeds and all seeds have the same period, and that when the period is maximal that there is just one sequence which visits all states exactly once. 

The generator has been run on Mac OSX systems with both x86_64 and PPC architecture. It compiles without problem on 64 bit Linux systems with appropriately configured compilers, but is expected to work correctly also on 32 bit systems wherever emulated 64 bit unsigned integer is available. 

It has been successfully compiled on several recent vintages of the GNU compiler gcc, and also on version 11 of the Intel compiler icc. 

The implementation uses bit shifts to implement the modular arithmetic, nevertheless, it produces the same output on big-endian and small-endian machines. This has been specifically tested on ppc. The program outputs uint64_t integers between 0 and 2^61-1 and therefore 61 lower bits are usable. The floating point output is configured to use long double floating point numbers which can mean different things on different architectures: on x86 and with gcc it is the extended precision 80-bit format; on PPC it is likely the same as the native 64 bit double so that there are differences in 21-22nd decimal digits between x86 and ppc. Keep in mind, though, that 61 bits of precision is good to only 18  decimal.

  The current implementation is not thread-safe. The facilities to allocate and initialize instances of the generator are very recent and less well tested compared to the rest. I am hoping to clarify this in some future release.



THEORY

Study of this recursion on a Galois field GF[p] is a complicated subject [6], suffice to say that for maximum period the characteristic polynomial of the matrix should have as many as possible of its eigenvalues among the primitive roots of the root-N-extended field GF[Root[p,N]]. There are no practical algorithms known to assure this in cases other than p=2. Recursions based on matrices whose characteristic polynomial is a primitive trinomial in GF[2^D] were advocated by Niederreiter [3] and have made their way into mainstream [7] in the form of Mersenne twisters of various D, for example D=19937. This is entirely due to the ease of finding primitive trinomials whenever 2^D-1  is a prime, namely it is sufficient that the polynomial is irreducible. The drawback is that one or more eigenvalues gets extremely close to the unit circle, and this gets worse (!) with large D.

   The eigenvalues of the matrix used by MIXMAX are well-separated from unity for all N. The present implementation uses by default the dimension N of the matrix equal to 32, and this can be changed by hand without adverse consequences. In practice, the terrible downfall of all recursive and matrix generators is their bad properties of distribution in the space of dimension equal or larger than 2*N. For example, the LCG generators suffer already in 2 dimensions, while Mersenne twisters have nearly 100% correlation between bits with certain lag, which is terrible indeed as in principle failure should only be apparent upon examining all 2*D bits. "Tempering" somewhat mitigates this situation. This brings us to the strongest point of the MIXMAX generator for physics simulations - if it is desired to integrate by Monte-Carlo sampling from a very high-dimensional space it is always possible to increase N such that it is greater than the dimensionality of that space. Fortunately, in the current implementation the computational complexity is of order O(N), rather than O(N^2) of the original or the O(N ln(N) ) of the improved implementation (ref. [6]). This means that per random number generated, the cost does not anymore grow with N. 

   One additional new feature, compared to other multiple recursive generators is the new function to skip over some large number of steps k by using the precomputed and stored matrix (A^k mod p) and therefore to have the very valuable ability to parallelize. This is implemented so far only on x86_64, and is for now dependent on endianness.

LITERATURE
In addition to the two papers where the generator was introduced, there exists among other the following literature which is most closely related.  Martin Luscher found in [5] the improvement to RCARRY which was sufficient in practice to overcome the high correlations intrinsic to the generators based on GF[2]. Luscher's RANLUX accomplishes this by skipping over the sequence which moves the eigenvalues further away from unit circle at the cost of speed.

[3]  H.Niederreiter. 
      "Finite fields, pseudorandom numbers, and quasirandom points," 
        in : Finite fields, Coding theory, and Advance in Communications and Computing. 
        (G.L.Mullen and P.J.S.Shine, eds) pp. 375-394, Marcel Dekker, N.Y. 1993.

[4] N.Z.Akopov, E.M.Madunts,A.B.Nersesian, G.K.Savvidy and W.Greiner, 
     "Fast K system generator of pseudorandom numbers."
      Proceedings of the XXVIII International Symposium Ahrenshoop, pp.281-286 (Wendisch-Rientz, Germany, 1994)
      arXiv:hep-lat/9311025

[5] M. Luscher, Computer Physics Communications  79 (1994) 100
     F. James,    Computer Physics Communications 79 (1994) 111
     
[6] "K-system generator of pseudorandom numbers on Galois field,"
     G. G. Athanasiu, E. G. Floratos, G. K. Savvidy
     DEMO-HEP 97/03, arXiv:physics/9703024
     International Journal of Modern Physics C, Volume 8, Issue 03, pp. 555-565 (1997).

[7] M. Matsumoto and T. Nishimura, 
    "Mersenne Twister: A 623-dimensionally equidistributed uniform pseudorandom number generator", 
    ACM Trans. on Modeling and Computer Simulation Vol. 8, No. 1, January pp.3-30 (1998) DOI:10.1145/272991.272995

[8] F.James, 
     "Finally, a theory of random number generation." 
     Fifth International Workshop on Mathematical Methods, CERN, Geneva, October 2001. 


TIMELINE

July, 1986
the original papers were out

July, 1987
Akopov and Savvidy talk to Fred James at CERN

December, 1991
the original papers are published in JCP

1994
Luscher publishes his method of improving the RCARRY generator by means of skipping

December 2008
I, Konstantin, find a way to implement the matrix recursion without explicitly using the matrix,
the computational complexity becomes O(N)

November 24, 2012
the initial version 0.01 is released on hepforge.org

Email me, Konstantin Savvidy, ksavvidis @at@ gmail.com
