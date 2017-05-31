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

int rayTrace = 1 ;

VelModel mp,ms ;
void printModel( char *text, double *x ) 
{
/* 	printf("Model %s : %9.4f %8.4f %8.4f %8.4f\n",text,x[0],x[1]*100.0,x[2]*100.0,x[3]) ; */
	printf("Model %s : %9.5f %9.5f %9.5f %9.5f\n",text,x[0],x[1],x[2],x[3]) ;
}
void insertRow( double *x, double *a, int i, int m, int n, double weight )
{
	while(n--) { a[i] = weight * *x++ ; i += m ; }
}
void testFit(char type, double x, double z, double t, double dx, double dz)
{
	double tt,dtdx,dtdz ;
	static	int count ;
/*
	double *model ;
	if( type == 'P' ) model = vFNModelP ;
	else model = vFNModelS ;
*/
	tt = vFtimeFromXZ(type,x,z,&dtdx,&dtdz) ;
	if( (count % 400) == 0 )fprintf(stderr,"xz = %10.4f %10.4f tt= %10.4f %10.4f dx %10.4f %10.4f dz %10.4f %10.4f\n"
		,x,z,t,tt,dtdx,dx,dtdz,dz) ;
}
int  locate( Solution *sol, Phase *pp )
{
	#define MAXPHASES 60 
	double x[4],x0[4],b[MAXPHASES],a[4*MAXPHASES] ;
	double dx[4], dxtest[4] ;
	double tt,dlon,distance,residual,sumP,sumS,stdP,stdS ;
	double rayp,dtdx,dxdp,km2lat,km2lon ;
	double dz,ttz,dtdz,dtdla,dtdlo ;
	double dl,dtdlax,dtdlox,damp, weight ;
	double kme,kmn, lenx, azi ;
	VelModel *vm ;
	Phase *p ;
	Station *s ;
	int np,i,j,iter,ns ;
	if (sol->index != pp->index ) {
		rLog(2,"Index in locate does not match",NULL) ;
		return 0 ;
	}
	p = pp ;
	while( p->index == sol->index ) p++ ;
	np = p-pp ;
	i = np ;
	while( i-- ) {
		p=pp+i ;
		s= p->statP ;
		distance = sDistance(sol->lat,s->lat,sol->lon - s->lon) ;
		if(distance > vFXMax * 0.9 ) np = i ;
	}
	if( np > MAXPHASES ) rLog(1,"locate: more than %d phases", (void*) MAXPHASES ) ;
	x0[0] = 0.0 ;         /* origin time, sec */
	x0[1] = sol->lat ;    /* latitude, degrees */
	x0[2] = sol->lon ;	/* longitude, degrees */
	x0[3] = sol->depth;	/* depth, km */
	dz = 0.005 ;
	km2lat = 6391*M_PI/180.0 ;
	km2lon = km2lat * cos(sol->lat * M_PI/180.0 ) ;
	for ( iter = 0 ; iter < 35 ; iter++) {
	  sumP = 0.0 ; sumS = 0.0 ; ns = 0 ;
	  if(shLogLevel > 4 ) printf("rayTrace = %d\n", rayTrace) ;
	  for( i = 0 ; i < np ; i++ ) {
		p = pp+i ;
		s = p->statP ;
		dlon = x0[2] - s->lon ;
		if( p->type == 'P' ) vm = &mp ; else vm = &ms ;
		if( p->type == 'P' ) weight = 1.0 ; else weight = 0.8 ;
		distance = sDistance(x0[1],s->lat,dlon) ;
		if( rayTrace ) { 
		  tt = timeFromDist(vm,distance,x0[3],&rayp,&dtdx,&dxdp) + x0[0] ;
		  ttz = timeFromDist(vm,distance,x0[3]+dz,&rayp,&dtdx,&dxdp) +x0[0] ;
		  dtdz = (ttz-tt)/dz ;
		} else tt = vFtimeFromXZ(p->type,distance,x0[3],&dtdx,&dtdz) ; 
/*		testFit(p->type,distance,x0[3],tt,dtdx,dtdz) ;  */
		dtdla = dtdx * ( x0[1] - s->lat )*km2lat*km2lat /  distance  ;
		dtdlo = dtdx * ( x0[2] - s->lon )*km2lon*km2lon /  distance  ;
		dx[0] = 1.0 ; dx[1] = dtdla ; dx[2] = dtdlo ; dx[3] = dtdz ;
		insertRow(dx,a,i,np,4,weight) ;
		b[i] = weight * ( p->pTime - tt ) ;
		residual = tt - p->pTime ;
		if( p->type == 'P' )sumP += residual*residual ; else {sumS +=  residual*residual ; ns++ ; }
/*
		dl = 0.0001 ;
		distance = gDistance(x0[1]+d2,s->lat,dlon) ;
		ttz = timeFromDist(vm,distan2e,x0[3],&rayp,&dtdx,&dxdp) +x0[0] ;
		dtdlax = (ttz-tt)/dl ;
		distance = gDistance(x0[1],s->lat,dlon+dl) ;
		ttz = timeFromDist(vm,distance,x0[3],&rayp,&dtdx,&dxdp) +x0[0] ;
		dtdlox = (ttz-tt)/dl ;

		printf("%s %c %10.6f %10.6f %10.6f %10.6f %10.6f ",s->name,p->type,p->pTime,tt,residual,distance,x0[3]) ;
		printf("%8.4f  %8.4f %8.4f ",dtdlax,dtdlox,1.0/dtdx) ;
		printModel(" dx ",dx) ; */
	   }
	   golubC(a,x,b,np,4) ;
	   kmn = x[1]*km2lat ; kme = x[2]*km2lon ;
	   lenx = sqrt( kmn*kmn + kme*kme + x[3]*x[3] ) ;
	   azi = atan2(kmn,kme) * 180.0/M_PI ; if( azi < 0 ) azi += 180 ;
	   damp = 20.0/lenx ;
	   damp = 2.6/lenx ;
#define DAMP 1.0
	   if(damp > DAMP) damp = DAMP ;
	   for( j = 0 ; j < 4 ; j++) x0[j] += damp*x[j] ;
	   if(shLogLevel >  3 )  {
	  	 printModel(" x ",x) ;
		 printModel(" x0",x0) ;
	   }
	   stdP = sqrt(sumP/(np-ns)) ;
	   stdS = sqrt(sumS/ns) ;
	   if(shLogLevel > 3 ) 
	     printf("iter = %2d  stdP =%9.6f stdS =%9.6f azi=%5.0f lenx=%7.3f damp=%7.2f\n",iter,stdP,stdS,azi,lenx,damp) ;
	   if( lenx < 0.01 ) break ;
	}
	sol->nP = np - ns ;
	sol->nS = ns ;
	sol->nIter = iter ;
	sol->sumP = sumP ;
	sol->sumS = sumS ;
	sol->length = lenx ;
	sol->stdP = stdP ;
	sol->stdS = stdS ;
	return np ;
}

#ifdef TEST

char *pModel = "silp.vel" ;
char *sModel = "sils.vel" ;
char *phaseFile = "../geysir/phase.dat" ;
char *solFile = "../geysir/ctloc2" ;
int skipEvents ;
int nEvents = 1 ;

doit( int skip )
{
	Phase *phases, *ip ;
	Solution *location, *lp ;
	int nPhases,nLoc ,i,j ;
	long long index,i2 ;
	double dist,azi ;
	Station *sp ;
	VelModel jm ;
/*	LocateStatus status ; */
	if( rayTrace ) {
	  initVelModel(20,&jm ) ;
/*  	  readVelModel(pModel,&jm) ; mp = resampleVelModel(&jm,1.00,50) ;
	  readVelModel(sModel,&jm) ; ms = resampleVelModel(&jm,1.00,50) ; */
	  initVelModel(20,&mp) ; readVelModel(pModel,&mp) ;
	  initVelModel(20,&ms) ; readVelModel(sModel,&ms) ; 
	} else vFInit() ;
	nPhases = readPhases(phaseFile,&phases ) ;
	nLoc = readCtloc(solFile,&location) ;
	lp = location + skip ;
	lp->nIter = 1 ;
	while( nEvents--) {
	  index = lp->index ;
	  printf("index=%ld\n",index) ;
	  ip = phases ;
	  while( ip->index < index ) ip++ ;
	  j = locate( lp, ip) ;
	  while( j-- ) {
		sp = ip->statP ;
		dist = sDistance(sp->lat,lp->lat,sp->lon-lp->lon ) ;
		azi = azAzimuth(lp->lat,sp->lat,sp->lon-lp->lon ) ;
		printf("%ld %s %c %10.6f %10.6f dist =%9.4f azi =%8.2f\n",ip->index,sp->name,ip->type,ip->pTime,ip->weight,dist,azi ) ;
		ip++ ;
	  }
	  printf("%d phases\n",ip-phases ) ;
	  lp++ ;
	}
}
int main(int ac, char **av) {
	int cc,n ;
	extern char *optarg ;
	extern int optind ;
	feenableexcept(FE_INVALID) ; 
	shLogLevel = 2 ;
	while( EOF != ( cc = getopt(ac,av,"l:p:s:fn:X:"))) {
	    switch(cc) {
	    case 'l' : skipEvents = atoi(optarg) ; break ;
	    case 'p' : phaseFile = optarg ; break ;
	    case 's' : solFile = optarg ; break ;
	    case 'f' : rayTrace = 0 ; break ;
	    case 'X' : vFXMax = atof(optarg) ; break ;
	    case 'n' : nEvents = atoi(optarg) ; break ;
	}}
	doit(skipEvents) ;
	return 0 ;
}
#endif
