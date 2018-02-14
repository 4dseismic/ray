

/* inverse errorfunction: 
https://en.wikipedia.org/wiki/Error_function#Inverse_functions */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "ray.h"
#define MK 100
double inverf( double z, int n )
{
	double ck[MK] ;
	int k,m ;
	double f,sum, fk2,fk ;
	if(n >= MK ) abort() ;
	fk =  sqrt(M_PI)*z/2.0 ; 
	fk2 = M_PI*z*z/4.0 ;
	ck[0] = 1.0 ;
	sum =  fk ;
	for( k = 1 ; k <= n ; k++) {
		ck[k] = 0.0 ;
		for( m = 0 ; m < k ; m++) {
			ck[k] += ( ck[m]*ck[k-1-m] ) / ( ( m+1.0 ) * ( m+m+1.0 )) ;
		}
		fk *= fk2 ;
		sum += ck[k] * fk / ( k+k+1.0 ) ;
		
	}
/*	for( k = 0  ; k <=n ; k++) printf("%12.8f ",ck[k] ) ;
	printf("\n") ; */
	return sum ;
}
double grandom( double x0, double x1, double xcenter, double sdev )
/* return a random number in interval from x0 to x1, 
	that has gaussian probabity density and is
	centered on xcenter with standard deviation sdev */
{
	double d0, d1, drand, xrand ;
	d0 = erf((x0-xcenter)/sdev) ;
	d1 = erf((x1-xcenter)/sdev) ;
	drand = d0 + random() * ( d1 - d0 ) / RAND_MAX ;
	xrand = xcenter + inverf(drand,90)* sdev ; 
	return xrand ;
}
#ifdef TEST
void test(double z, int n )
{
	double sum ;
	sum = inverf(z,n) ;
	printf("%3d %12.8f %12.8f %12.8f\n",n,z,sum,erf(sum) ) ;
}
void testGrandom()
{
	int i ;
	double xr ;
	for( i = 0 ; i < 10000 ; i++) {
		xr = grandom( 3.0,5.0, 4.5, 0.4 ) ;
		printf("%8.4f\n",xr) ;
	}
}
void testInverf() 
{
	double e ;
	test(0.5,5) ;
	test(0.5,6) ;
	test(0.5,7) ;
	test(0.5,8) ;
	test(0.9,5) ;
	test(0.9,6) ;
	test(0.9,7) ;
	test(0.9,8) ;
	test(0.9,20) ;
	test(-0.9,20) ;
	test(0.99,20 ) ;
	test(0.99,99 ) ;
}
int main(int ac , char **av) 
{
/*	testInverf() ; */
	testGrandom() ;
	return 0 ;
}
#endif 
