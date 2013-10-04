#include <stdio.h>
#include <math.h>
#include "mixmax2008.h"

int main( int argc, /* number of arguments */ char *argv[] /* arguments */){
    int j,k;
	unsigned int p;

	myuint *X=rng_alloc();
    //myuint X[N+1];

	fprintf(stderr,"Welcome to the Yerevan random number generator!\nThe curent matrix size is %u\n"
		   "(the actual matrix is not kept in memory in this new efficient implementation)",N);
    fprintf(stderr,"\nWe will now generate a power of our matrix A^m\n "); 
    fprintf(stderr,"\nHow many iterations m do you want?\n(2^10=1024, 2^20=1048576, 2^30=1073741824)\nEnter m: "); scanf("%u", &p);

		printf("\nHere is convenient Mathematica Input:\n AT={");
		for (k=1; (k<=N); k++) {
			seedvector(X,k);
			for (j=0; (j<p); j++) {
				IT(X);
			}
			printf("\n{");
			printf("%llu", X[1] );
			for (j=2; (j<=N); j++) {
				printf(",\t%llu", X[j] );
			}
			printf("},");
		}
		printf("\b} \n A = Transpose[AT] \n d = Det[A] \n AM = Mod[A, 2^61-1] \n cp = CharacteristicPolynomial[AM, x] \n "
			   "N[Solve[cp == 0, x]] \n N[Eigenvalues[AM]] \n Abs[%%]\n"
				"M1K = Mod[MatrixPower[A, 2^10], 2^61 - 1]\n "); 
		//\n Maxima input: eq:charpoly(N,x);\n p:radcan(eq); r:allroots(p);
	printf("\n\nok\n");
	rng_free(X);
	return 1;
}
