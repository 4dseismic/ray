
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pr.h"
#define K (RAD2DEG*60.0)

typedef struct { double l,b ; } LatLon ;
typedef struct { char map[32] ; double lon,lat,y,x ; } MapPoint ;

V3 prPCenter, prPN1,prPN2 ;
double prPixPro[6], prIPixPro[6]  ;

int prMapScale=25 ;

void prLatLon2V3( V3 *v, LatLon *l )
/* convert from lat-long (degrees) to direction vector */
{
	double co ;
	v->z = sin(l->b) ;
	co = cos(l->b) ;
	v->x = cos(l->l)*co ;
	v->y = sin(l->l)*co ;
}
void prMakeV3(V3 *c, double lon, double lat )
/* make a vector from lat - lon values */
{
	LatLon a ;
	a.l = lon ;
	a.b = lat ;
	prLatLon2V3(c,&a) ;
}
double prDotProd( V3 *a, V3 *b)
/* vector dot or scalar product */
{
	return a->x * b->x + a->y * b->y + a->z * b->z ;
}
void prXProd( V3 *c, V3 *a, V3 *b)
/* outer or cross prodoct of vectors a and b */
{
	c->x = a->y * b->z - a->z * b->y ;
	c->y = a->z * b->x - a->x * b->z ;
	c->z = a->x * b->y - a->y * b->x ;
}
double prV3Length( V3 *a)
{
	return sqrt(a->x * a->x  +  a->y * a->y  +  a->z * a->z) ;
}
void prV3Normalize( V3 *b, V3 *a)
/* place vector a normalized in b, may overlap */
{
	double rr ;
	rr =  1.0/sqrt(a->x * a->x  +  a->y * a->y  +  a->z * a->z) ;
	b->x = rr*a->x ;
	b->y = rr*a->y ;
	b->z = rr*a->z ;
}
void prScalMul2( V3 *c, double co, V3 *a, double si, V3 *b)
{
	c->x = co * a->x + si * b->x ;
	c->y = co * a->y + si * b->y ;
	c->z = co * a->z + si * b->z ;
}
double prAngle( V3 *a, V3 *b)
/* angle between two unitvectors, range 0 to 180 */
{	
	V3 n ;
	double co,si ;
	co = prDotProd(a,b) ;
	prXProd(&n,a,b) ;
	si = prV3Length(&n) ;
	return atan2(si,co);
}
double prAzimuth( V3 *a, V3 *b) 
/* compute azimuth from a to b, <-180-180 degrees> */
{
	V3 n1,n2,n ;
	double angle,co,si ;
	n1.x = a->y ;
	n1.y = -a->x ;
	n1.z = 0.0 ;
	prXProd(&n2,a,b) ;
	co = prDotProd(&n1,&n2) ;
	prXProd(&n,&n2,&n1) ;
	si = prDotProd(&n,a) ;
	angle = atan2(si,co) ;
	return angle ;
}
void prProject( double *x, double *y, double lon, double lat)
/* project from lat-lon to cartesian using prjection center defined by
 *  prSetProjCenter, result is in arc-minutes or nautical miles. */
{
	V3 b,n2 ;
	double r,si,co ;
	prMakeV3(&b,lon,lat) ;
	r = prAngle(&b,&prPCenter) ;
	prXProd(&n2,&b,&prPCenter) ;
	co = prDotProd(&prPN1,&n2) ;
	si = prDotProd(&prPN2,&n2) ;
	if( r > 0.0 ) r = K*r / sqrt(si*si + co*co ) ;
	*x = -si*r ;
	*y = co*r ;
/*	printf("x=%12.6f y=%12.6f\n",*x,*y) ; */
}
V3 prDest(V3 *a, double dist, double azimuth )
/* return point in direction azimuth and dist away from a */
{
	double cosd,sind,cosa,sina ;
	V3 n1,n2,n3,res ;
	n1.x = -(a->y) ;
	n1.y = (a->x) ;
	n1.z = 0.0 ;
	prV3Normalize(&n1,&n1) ;
	prXProd(&n2,a,&n1) ;
	cosd = cos( dist ) ;
	sind = sin( dist ) ;
	cosa = cos( azimuth ) ;
	sina = sin( azimuth ) ;
	prScalMul2(&n3,cosa,&n2,sina,&n1) ;
	prScalMul2(&res,cosd,a,sind,&n3) ;
#ifdef TEST
	printf("dist=%g azimuth=%g\n",dist,azimuth) ;
	printf("a  = %s %g\n",prLatLon2Ascii(a),prV3Length(a)) ;
	printf("n1 = %s %g\n",prLatLon2Ascii(&n1),prV3Length(&n1)) ;
	printf("n2 = %s %g\n",prLatLon2Ascii(&n2),prV3Length(&n2)) ;
	printf("n3 = %s %g\n",prLatLon2Ascii(&n3),prV3Length(&n3)) ;
	printf("res  = %s %g\n",prLatLon2Ascii(&res),prV3Length(&res)) ;
#endif
	return res ;
}
void prProjectIP( double *lon, double *lat, double xp, double yp)
/* project from pixels to lat-lon , using projection set by prSetProjCenter
 * and pixel projection defined by prIFitPoint() */
{
	double x,y,r,rr,a ;
	V3 v ;
	x = prIPixPro[0]  +  xp * prIPixPro[1]  +  yp * prIPixPro[2] ;
	y = prIPixPro[3]  +  xp * prIPixPro[4]  +  yp * prIPixPro[5] ;
	r = sqrt(x*x + y*y) / K ;
	a = atan2(x,y) ;
	v = prDest(&prPCenter,r,a) ;
	*lon = atan2(v.y,v.x) ;
	rr = sqrt( v.x * v.x + v.y * v.y ) ;
	*lat = atan2( v.z,rr) ;
}
void prProjectP( double *xp, double *yp, double lon, double lat)
/* project from lat-lon to pixels, using projection set by prSetProjCenter
 * and pixel projection defined by prFitPoint() */
{
	double x,y ;
	prProject(&x,&y,lon,lat) ;
	*xp = prPixPro[0]  +  x * prPixPro[1]  +  y * prPixPro[2] ;
	*yp = prPixPro[3]  +  x * prPixPro[4]  +  y * prPixPro[5] ;
}
void prCheckFit( MapPoint *m, int n) 
{
	double l,b,x,y,xp,yp ;
	int i ;
	for( i = 0 ; i < 6 ; i++) printf("%g ",prPixPro[i]) ;
	printf("\n") ;
	for( i = 0 ; i < n ; i++) {
		prProjectP(&xp,&yp,m->lon,m->lat) ;
		printf("%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f\n",
			RAD2DEG*m->lon,RAD2DEG*m->lat,m->x,m->y,xp,yp) ;
		m++ ;
	}
}
void prIsbladMap(  double lon, double lat, char *buffer )
{
	double l[]= {-22.9, -22.7417, -22.6033, -19.0179, -19.0187,
		-19.025, -15.139, -15.2958, -15.4392 } ;
	double b[] = { 66.0167, 64.9458, 63.8683, 66.0638, 64.9888,
		63.9167, 66.0183, 64.9417, 63.8675 } ;
	double dm,d,ll,bb ;
	int i,im ;
	ll = lon*RAD2DEG ;
	bb = lat*RAD2DEG ;
	dm =10000.0 ;
	for( i = 0 ; i < 9 ; i++ ) {
		d = fabs(ll - l[i]) + fabs(bb - b[i]) ;
/*		printf("i=%d d=%g  place=%s\n",i,d,nmLatLon2Ascii(lon,lat)) ;*/
		if( d <= dm ) {
			im = i ; 
			dm = d ;
		}
	}
	sprintf(buffer,"isblad%d",im+1) ;
}
#define N 6 
#define M 8
void prInvFitPoint( MapPoint *m, int n )
/* fit transformation from pixels to projected coordinates*/
{
	double a[N*M], b[M] ;
	int i,j ;
	double xx,yy ;
	MapPoint *p ;
	memset(a,0,N*M*sizeof(double)) ;
	memset(b,0,M*sizeof(double)) ;
	j = 0 ;
	for( i = 0 ; i < n ; i++ ) {
		p = m+i ;
		prProject(&xx,&yy,p->lon,p->lat) ;
		b[j] = xx ;
		a[j     ] = 1.0 ;
		a[j +  M] = p->x ;
		a[j +2*M] = p->y ;
		j++ ;
		b[j] = yy ;
		a[j +3*M] = 1.0 ;
		a[j +4*M] = p->x ;
		a[j +5*M] = p->y ;
		j++ ;
	}
	golubC(a,prIPixPro,b,M,N) ;
}
#ifdef TEST
char prTestDmaMap()
{
	prDmaMap( DEG2RAD * -19.9, DEG2RAD * 64.1 ) ;
	prDmaMap( DEG2RAD * -19.3, DEG2RAD * 64.1 ) ;
	prDmaMap( DEG2RAD * -19.3, DEG2RAD * 64.3 ) ;
	prDmaMap( DEG2RAD * -19.9, DEG2RAD * 64.3 ) ;
	prDmaMap( DEG2RAD * -19.9, DEG2RAD * 64.6 ) ;
	prDmaMap( DEG2RAD * -19.3, DEG2RAD * 64.6 ) ;
}
int main(int ac, char **av)
{
	V3 center, car, dest ;
	double d,a,l,b,ll,bb ;
	double x,y ;
	l = DEG2RAD*atof(av[1]) ;
	b = DEG2RAD*atof(av[2]) ;
	prMakeV3(&center,-22*DEG2RAD,64*DEG2RAD) ;
	prMakeV3(&car,l,b) ;
	d = prAngle(&center,&car) ;
	a = prAzimuth(&center,&car) ;
	prSetProjCenter(-22*DEG2RAD,64*DEG2RAD) ;
	prProject(&x,&y,l,b) ;
	dest = prDest(&center,d,a) ;
	printf("%12.4f %12.4f %12.4f %12.4f %s\n",l*RAD2DEG,b*RAD2DEG,
			d*RAD2DEG*60,a*RAD2DEG,prLatLon2Ascii(&dest)) ;
	printf("x=%12.5f y=%12.5f %s \n",60*RAD2DEG*x,60*RAD2DEG*y,
		prLatLon2Ascii(&prPN2)) ;
	prReadMapList("point.list","is1613_3") ;
	prProjectP(&x,&y,l,b) ;
	printf("%g %g, %s\n",x,y,nmLatLon2Ascii(l,b)) ;
	prProjectIP(&ll,&bb,x,y) ;
	printf("%g %g, %s\n",x,y,nmLatLon2Ascii(ll,bb)) ;
	prTestDmaMap() ;
	return 0 ;
}
#endif

