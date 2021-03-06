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

/* #define MAXPAR 45 */
extern VelModel mp,ms ;
double vFNModelP[MAXPAR],  vFNModelS[MAXPAR] ;
int vFNPar, vFOrder ;
typedef enum { poly, ft} FitType ;
FitType vFType ;
double vFdfx , vFdfz ;

double vFXMax = 100.0 ;
double vFZMax = 15.0  ;
int nX = 50 ;
int nZ = 20 ;

double dfreqz, dfreqx 	;

typedef struct { double x,z,t,tfit ; } TimePoint ;
TimePoint *vList ;
int nVList ;

char *velFile = "silp.vel" ;

void vFPrintModel( char *text, double *x, int n ) 
{
	int i ;
	printf("Model %s vFNPar=%d vFOrder=%d\n",text,vFNPar,vFOrder ) ;
	for( i = 0 ; i < n ; i++ ) printf("%9.4f ",x[i]) ;
	printf("\n") ;
}
void printModelE( char *text, double *x, int n ) 
{
	int i ;
	printf("Model %s vFNPar=%d vFOrder=%d\n",text,vFNPar,vFOrder ) ;
	for( i = 0 ; i < n ; i++ ) printf("%20.15e ",x[i]) ;
	printf("\n") ;
}
void vFInsertRow( double *x, double *a, int i, int m, int n, double weight )
{
	while(n--) { a[i] = weight * *x++ ; i += m ; }
}
void travelTTable( char * vfile)
{
	VelModel vm,jm ;
	int ix,iz ;
	double x,z,dx,dz,t,dtdx,p,dxdp ;
	TimePoint *tp ;
	initVelModel(20,&jm) ;
	nVList = nX*nZ ;
	if( NULL == vList ) vList = ( TimePoint *) malloc( nVList * sizeof(TimePoint) ) ;
	tp = vList ;
	readVelModel(vfile,&jm) ;
	vm = resampleVelModel(&jm,0.25,50) ;
	dz = vFZMax/nZ ; dx = vFXMax/nX ;
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
void travelTTableNoIter( char * vfile)
{
	VelModel vm,jm ;
	int ix,iz ;
	double x,z,dx,dz,t,dtdx,p,dxdp ;
	TimePoint *tp ;
	initVelModel(20,&jm) ;
	nVList = nX*nZ ;
	if( NULL == vList ) vList = ( TimePoint *) malloc( nVList * sizeof(TimePoint) ) ;
	tp = vList ;
	readVelModel(vfile,&jm) ;
	vm = resampleVelModel(&jm,0.25,50) ;
	dz = vFZMax/nZ ; dx = vFXMax/nX ;
        for( iz = 0 ; iz < nZ ; iz++) {
	   z = (0.5 + iz) * dz ;
	   for( ix = 0 ; ix < nX ; ix++ ) {
	 	x = (0.5 + ix) * dx ;
		tp->x = x ; tp->z = z ;
		tp->t = timeFromDist(&vm,x,z,&p,&dtdx,&dxdp) ;
		tp++ ;
	    }
	}
}
void makeTTable( VelModel *vmp )
{
	VelModel vm ;
	int ix,iz ;
	double x,z,dx,dz,t,dtdx,p,dxdp ;
	TimePoint *tp ;
	printf("Enter makeTTable \n") ;
	if(shLogLevel > 5 )printf("entering makeTTable, v[0]= %10.4f\n",vmp->v[0] ) ;
	nVList = nX*nZ ;
	if( NULL == vList ) vList = ( TimePoint *) malloc( nVList * sizeof(TimePoint) ) ;
	tp = vList ;
	vm = resampleVelModel(vmp,0.25,50) ;
	dz = vFZMax/nZ ; dx = vFXMax/nX ;
	for( ix = 0 ; ix < nX ; ix++ ) {
	    x = (0.5 + ix) * dx ;
	    for( iz = 0 ; iz < nZ ; iz++) {
		z = (0.5 + iz) * dz ;
		tp->x = x ; tp->z = z ;
		tp->t = timeFromDist(&vm,x,z,&p,&dtdx,&dxdp) ;
		tp++ ;
	    }
	}
	printf("Leave makeTTable \n") ;
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
	d[7] = 1.0/(x-z + vFZMax+4.0) ;
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
/*
double xfac = 0.02 ;
double zfac = 0.20 ;
*/
double xfac = 1.0 ;
double zfac = 1.0 ;
double timeFuncP( double x, double z, double *c, double *d ) 
{
/* 2d polonomial  nPar = order * (order+1)/2  */
	double sum,xx,zz, o ;
	int j,i,ij ;
	ij = 0 ;
	zz = 1.0 ; 
	x *= xfac ;
	z *= xfac ;
	for( i = 0 ; i < vFOrder ; i++ ) {
		xx = 1.0 ;
		for( j = 0 ; j < vFOrder-i ; j++ ) {
			d[ij++] = xx*zz ;
			xx *= x ;
		}
		zz *= z ;
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
double vFtimeFromXZ(char type, double x, double z, double *dtdx, double *dtdz)
{
	double sum, sumdx, sumdz ;
	double xx,zz,xx0,zz0 ;
	int j,i,ij ;
	double f ;
	double *model ;
	if( type == 'P' ) model = vFNModelP ;
	else model = vFNModelS ;
	sum = 0.0 ; sumdx = 0.0 ; sumdz = 0.0 ;
	ij = 0 ;
	zz = 1.0 ; zz0 = 1.0 ;
	x *= xfac ;
	z *= zfac ;
	for( i = 0 ; i < vFOrder ; i++ ) {
		xx = 1.0  ;  xx0 = 1.0 ; 
		for( j = 0 ; j < vFOrder-i ; j++ ) {
			f = model[ij++] ;
			sum += f*xx*zz ;
			sumdx += j*f*xx0*zz ;
			sumdz += i*f*xx*zz0 ;
			xx0 = xx ;
			xx *= x ;
		}
		zz0 = zz ;
		zz *= z ;
	}
	*dtdx = sumdx*xfac ;
	*dtdz = sumdz*zfac ;
	return sum ;
}
void linearFit(double *c)
{
	int i,m ;
	static double *a, *b ;
	double d[MAXPAR],dc[MAXPAR] ;
	double t,w,sum ;
	TimePoint *tp ;
	printf("Entering linearFit nVList = %d\n",nVList ) ;
/*	timeFuncSet(poly,29,0.0,0.0) ; */
	timeFuncSet(poly,29,130.0,20.0) ;
	m = nVList ;
	if( NULL == a ) {
		a = calloc(MAXPAR*m,sizeof(double)) ;
		b = calloc(       m,sizeof(double)) ;
	}
	tp = vList ;
	w = 1.0 ;
	for( i = 0 ; i < m ; i++) {
		t = timeFunc( tp->x, tp->z,c,d) ;
		vFInsertRow(d,a,i,m,vFNPar,w) ;
		b[i] = w*tp->t ;
		tp++ ;
	}
	golubC(a,c,b,m,vFNPar) ;
	if(shLogLevel > 5 ) vFPrintModel("c",c,vFNPar) ;
	tp = vList ;
	sum = 0 ; 
	for( i = 0 ; i < m ; i++) {
		tp->tfit = timeFunc( tp->x, tp->z,c,d) ;
		t = tp->t - tp->tfit ;
		sum += t*t ;
		tp++ ;
	}
	if( shLogLevel > 4 )printf("rms deviation: %8.3f\n", sqrt(sum/m) ) ;
}
void vFInitFromMemory()
{	/* velocity function is already in memory */
	makeTTable(&mp) ;
	linearFit(vFNModelP ) ;
	makeTTable(&ms) ;
	linearFit(vFNModelS ) ;
}
void vFInit()
{	/* velocity function is read from disk */
	int iterate;
	iterate = 0 ;
	if (iterate ) travelTTable("silp.vel") ;
	else travelTTableNoIter("silp.vel") ;
	linearFit(vFNModelP) ;
	if ( iterate )travelTTable("sils.vel") ;
	else travelTTableNoIter("sils.vel") ;
	linearFit(vFNModelS) ;
}
#ifdef TEST
VelModel mp,ms ;
int main(int ac, char **av) {
	int cc,n ;
	feenableexcept(FE_INVALID) ; 
	shLogLevel = 2 ;
	while( EOF != ( cc = getopt(ac,av,"f"))) {
	    switch(cc) {
	}}
	vFInit() ;
	dumpTimeTable() ;
	return 0 ;
}
#endif
