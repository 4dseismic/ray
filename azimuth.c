#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pr.h"
#include "ray.h"

double deg2rad = M_PI/180.0 ;

void printV( V3 *a)
{
	double l,b,z ;
	z = sqrt( a->x * a->x + a->y * a->y ) ;
	b = atan2(a->z,z) ;
	l = atan2(a->y,a->x) ;
	printf("(%10.6f %10.6f %10.6f) %11.6f  %11.6f\n",
		a->x,a->y,a->z,l/deg2rad,b/deg2rad) ;
}

double azAzimuth( double la1, double la2, double dlon )
/* return azimuth, all values are in degrees */
{
	V3 a,b ;
	double b1,b2,dlo,azim ;
	b1 = la1 * deg2rad ;
	b2 = la2 * deg2rad ;
	dlo = dlon * deg2rad ;
	prMakeV3(&a,0.0,b1) ;
	prMakeV3(&b,dlo,b2) ;
/*	printV(&a) ; printV(&b) ; */
	azim = prAzimuth(&a,&b) ; 
	return azim/deg2rad ;
}
#ifdef TEST
main()
{
	double jj ;
	jj = azAzimuth( 64,63,0.0) ;
	printf("jj=%10.4f\n") ;
}
#endif
