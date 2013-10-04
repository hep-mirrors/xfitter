#define N 32 //128
/* For N=32, will generate (and print) a billion numbers in under 4min, on Intel Xeon 2.8Ghz . 
 Printing takes 75% of that time. 
 */
#define BITS 61

#include <stdio.h>
#include "mixmax2008.h"


int main( int argc, /* number of arguments */ char *argv[] /* arguments */){
    int j,k;
	myuint *X=rng_alloc();
  
	fprintf(stderr,"Welcome to the Yerevan random number generator!\nThe curent matrix size is %u\n"
		   "(the actual matrix is not kept in memory in this new efficient implementation)",N);

		printf("\nHere are the unit vectors after acting by 1M skip:\n {");
		for (k=1; (k<=N); k++) {
			seedvector(X,k);
			skip_1M(X);
			printf("\n{");
			printf("%llu", X[1] );
			for (j=2; (j<=N); j++) {
				printf(",\t%llu", X[j] );
			}
			printf("},");
		}
		printf("\b}\n"); 

		
		printf("\nHHere are the unit vectors after acting by 1G skip:\n {");
		for (k=1; (k<=N); k++) {
			seedvector(X,k);
			skip_1G(X);
			printf("\n{");
			printf("%llu", X[1] );
			for (j=2; (j<=N); j++) {
				printf(",\t%llu", X[j] );
			}
			printf("},");
		}
		printf("\b}\n"); 
	return 1;
}
