
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "math.h"
int shuffle( int n ) 	/* return number from a shufled list of integers [0..n-1] */
{
	static int nn, nIn, *list ;
	int i,j,k,t ;
	if ( nn != n ) {
		if( list ) free(list) ;
		list = malloc(n*sizeof(int)) ;
		for( i = 0 ; i < n ; i++ ) list[i] = i ;
		nIn = 0 ; 
		nn = n ;
	}
	if(nIn == 0 ) {	/* shuffle when needed */
		i = 2*n ;
		while( i-- ){
			k = random() % n ;	
			j = random() % n ;
			t = list[k] ;
			list[k] = list[j] ;
			list[j] = t ;
		}	
		nIn = n ;
	}
	nIn-- ;
	return list[nIn] ;
}
#ifdef TEST 
int main(int ac, char **av) 
{
	int i,n ;
	time_t tt  ;
	srandom(time(&tt)) ;
	n = 6 ;
	for( i = 0 ; i < 100 ; i++) {
		printf(" %2d | ",shuffle(n) ) ;
		if( i%n == n-1 ) printf("\n" ) ;
	}
	return 1 ;
}
#endif
