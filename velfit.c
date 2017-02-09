/*
	Routines to solve for eq locations.
	eik dec 2016
*/
#include <stdlib.h>
#include <stdio.h>
#define _GNU_SOURCE
#include <fenv.h>
#include <math.h>
#include "ray.h"

int vFNPar, vFOrder ;
typedef enum { poly, ft} FitType ;
FitType vFType ;
double vFdfx , vFdfz ;

double xMax = 100.0 ;
double zMax = 15.0  ;
int nX = 50 ;
int nZ = 20 ;

double dfreqz, dfreqx ;

typedef struct { double x,z,t,tfit ; } TimePoint ;
TimePoint *vList ;
int nVList ;

char *velFile = "silp.vel" ;

void printModel( char *text, double *x, int n ) 
{
	int i ;
	printf("Model %s ",text ) ;
	for( i = 0 ; i < n ; i++ ) printf("%9.4f ",x[i]) ;
	printf("\n") ;
}
void insertRow( double *x, double *a, int i, int m, int n, double weight )
{
	while(n--) { a[i] = weight * *x++ ; i += m ; }
}
void travelTTable()
{
	VelModel vm,jm ;
	int ix,iz ;
	double x,z,dx,dz,t,dtdx,p,dxdp ;
	TimePoint *tp ;
	initVelModel(20,&jm) ;
	nVList = nX*nZ ;
	vList = ( TimePoint *) malloc( nVList * sizeof(TimePoint) ) ;
	tp = vList ;
	readVelModel(velFile,&jm) ;
	vm = resampleVelModel(&jm,0.25,50) ;
	dz = zMax/nZ ; dx = xMax/nX ;
	for( ix = 0 ; ix < nX ; ix++ ) {
	    x = (0.5 + ix) * dx ;
	    for( iz = 0 ; iz < nZ ; iz++) {
		z = (0.5 + iz) * dz ;
		tp->x = x ; tp->z = z ;
		tp->t = timeFromDist(&vm,x,z,&p,&dtdx,&dxdp) ;
		tp++ ;
	    }
	}
}
void dumpTimeTable()
{
	int i,j ;
	TimePoint *tp ;
	FILE *of ;
	tp = vList ;
	of = fopen("time.dump","w") ;
	for( i = 0 ; i < nX ; i++) {
		for( j = 0 ; j < nZ ; j++) {
			fprintf(of,"%7.3f",tp->t - tp->tfit) ;
			tp++ ;
		}
		fprintf(of,"%3d %5.1f\n",i,(tp-1)->x) ;
	}
	fclose(of) ;
}
double timeFunc4( double x, double z, double *c, double *d ) 
{
	double sum ;
	int j ;
	d[0] = 1.0 ;
	d[1] = x ;
	d[2] = z ;
	d[3] = x*z ;
	d[4] = 1.0/(x+2.5) ;
	d[5] = 1.0/(z+3.0) ;
	d[6] = 1.0/(x+z+3.0) ;
	d[7] = 1.0/(x-z+zMax+4.0) ;
	sum = 0.0 ;
	for( j = 0 ; j < vFNPar ; j++) { sum += c[j] * d[j] ; }
	return sum ;
}
double timeFuncSet( FitType type , int nPar, double xSpan,double zSpan ) 
/* use Fourier cosine transform to parameterize traveltime */
{
	double o, dfx, dfz ;
	int order ;
	vFType = type ;
	vFNPar = nPar ;
	
/*
	nz = order ;
	nx = 2 * order ; */
	if( type == poly ) {
		o = ( sqrt(8*nPar + 1.0 ) - 1.0 )*0.5 ;
		vFOrder = ceil(o) ;
		return ;
	}
	if( type == ft ) {
		vFOrder = ceil( sqrt( 1.0 * nPar ) ) ;
		vFdfx = 2.0 * M_PI / ( 3.0*xSpan ) ;
		vFdfz = 2.0 * M_PI / ( 3.0*zSpan ) ;
		return ;
	}
}
double timeFuncF( double x, double z, double *c, double *d ) 
{
	int i,j,ij ;
	double afx,afz, sum ;
	double cx, cx0, cx1, cz, cz0,cz1 ;
	ij = 0 ;
	afx = 2.0 * cos(x* vFdfx ) ;
	afz = 2.0 * cos(z* vFdfz ) ;
	cx0 = 1.0 ; cz0 = 1.0 ; cx1 = 1.0 ; cz1 = 1.0 ;
	for( i = 0 ; i < vFOrder ; i++ ) {
		cz = afz * cz1 - cz0 ;
		for ( j = 0 ; j < vFOrder ; j++ ) {
			cx = afx * cx1 - cx0 ;
			d[ij++] = cz*cx ;
			cx0 = cz1 ;
			cx1 = cx ;
		}	
		cz0 = cz1 ;
		cz1 = cz ;
	}
	sum = 0.0 ;
	for( j = 0 ; j < vFNPar ; j++) { sum += c[j] * d[j] ; }
	return sum ;
}
double timeFuncP( double x, double z, double *c, double *d ) 
{
/* 2d polonomial  nPar = order * (order+1)/2  */
	double sum,xx,zz, o ;
	int j,i,ij ;
	ij = 0 ;
	zz = 1.0 ; 
	for( i = 0 ; i < vFOrder ; i++ ) {
		xx = 1.0 ;
		for( j = 0 ; j < vFOrder-i ; j++ ) {
			d[ij++] = xx*zz ;
			xx *= x*0.02 ;
		}
		zz *= z*0.20 ;
	}
	sum = 0.0 ;
	for( j = 0 ; j < vFNPar ; j++) { sum += c[j] * d[j] ; }
	return sum ;
}
double timeFunc(double x, double z, double *c, double *d ) 
{
	if( vFType == poly ) return timeFuncP( x, z, c, d ) ;
	if( vFType == ft ) return timeFuncF( x, z, c, d ) ;
/*	return timeFunc4( x, z, c, d ) ; */
}
void linearFit()
#define MAXPAR 45
{
	int i,m ;
	double *a, *b,c[MAXPAR],d[MAXPAR],dc[MAXPAR] ;
	double t,w,sum ;
	TimePoint *tp ;
/*	timeFuncSet(poly,29,0.0,0.0) ; */
	timeFuncSet(poly,29,130.0,20.0) ;
	m = nVList ;
	a = calloc(MAXPAR*m,sizeof(double)) ;
	b = calloc(       m,sizeof(double)) ;
	tp = vList ;
	w = 1.0 ;
	for( i = 0 ; i < m ; i++) {
		t = timeFunc( tp->x, tp->z,c,d) ;
		insertRow(d,a,i,m,vFNPar,w) ;
		b[i] = w*tp->t ;
		tp++ ;
	}
	golubC(a,c,b,m,vFNPar) ;
	printModel("c",c,vFNPar) ;
	tp = vList ;
	sum = 0 ; 
	for( i = 0 ; i < m ; i++) {
		tp->tfit = timeFunc( tp->x, tp->z, c,d) ;
		t = tp->t - tp->tfit ;
		sum += t*t ;
		tp++ ;
	}
	printf("rms deviation: %8.3f\n", sqrt(sum/m) ) ;
}
int main(int ac, char **av) {
	int cc,n ;
	feenableexcept(FE_INVALID) ; 
	shLogLevel = 2 ;
	while( EOF != ( cc = getopt(ac,av,"f"))) {
	    switch(cc) {
	}}
	travelTTable() ;
	linearFit() ;
	dumpTimeTable() ;
	return 0 ;
}
