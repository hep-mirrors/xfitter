#include <stdio.h>
#include "mixmax2008.h"

int main( int argc,  char *argv[] ){
	myuint j,p;
	fprintf(stderr,"Welcome to the Yerevan random number generator!\nThe curent matrix size is %u\n"
		   "(the actual matrix is not kept in memory in this new efficient implementation)",N);

	fprintf(stderr,"\nHow many numbers?\n(2^10=1024, 2^20=1048576, 2^30=1073741824)\nEnter m: "); scanf("%llu", &p);
	myuint *S=rng_alloc();
	seedvector(S,1);   // seed with unit vector 
	//seed_lcg(S);     // or try seeding with some arbitrary lcg

	for (j=1; j<=p ; j++) {
		//z=get_next_float(S);
		//printf("%1.20LF\n", z );             // for floating point number on [0,1)
		printf("  %19llu\n", get_next(S) );  // for integer in 0 .. 2^61-2
		//printf("  %19llu\n", j );              // for timing, dont compute anything, just print the counter j
	}
		
	
	printf("ok\n");
	rng_free(S);
	return 1;
}
