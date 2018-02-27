
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>
#define _GNU_SOURCE
#include <fenv.h>
#include "../rayt/ray.h"

Event *events ;
Phase *phases ;
Solution *solutions ;
int nEvent, nPhase, nSol, nIndexList ;
extern VelModel mp,ms ;
VelModel initP,initS ;
long long *indexList ;
int gaussFlag, metropolis,noRandomize ;

/*		Default values */
int nIterations = 500 ;
char *velPrefix = "sil" ;
char *baseName = "testing" ;
int phasesMin = 8 ;
int phasesMax = 999 ;
int minSorP = 2 ;  /* minimum number of either p or s phases */
double zMax = 15.0 ;
double zMin = 4.0 ; 
double latMin = 63.0 ;
double latMax = 65.0 ;
double lonMin = -22.0 ;
double lonMax = -19.4 ; 
double distMax  = 90.0 ;
long long indexMin = 19910700000000000 ;
long long indexMax = INDEXEND ;
double weightS = 0.366 ;  
double temperature  =  0.03 ;
char *indexFile ;
int indexFileFraction = 100 ; /* precentage of entries from indexFile to use */
int nLayers = 7 ;
int skipLayers = 0 ;

void printParameters( FILE *f)
{
	fprintf(f,"velPrefix = %s baseName = %s\n",velPrefix,baseName) ;
	fprintf(f,"phases from %d to %d minimum S or P phases %d\n",
		phasesMin,phasesMax,minSorP) ;
	fprintf(f,"depths from %8.3f to %8.3f\n",zMin,zMax) ;
	fprintf(f,"latitude from %8.4f to %8.4f\n",latMin,latMax) ;
	fprintf(f,"longitude from %8.4f to %8.4f\n",lonMin,lonMax) ;
	fprintf(f,"index from %ld to %ld\n",indexMin,indexMax) ;
	fprintf(f,"nLayers = %d, skipLayers=%d, distMax=%8.2f\n",nLayers,skipLayers,distMax) ;
}
void printSolutions() 
{
	int i ;
	Solution *p ;
	for( i = 0 ; i < nSol ; i++) {
	  p = solutions + i ;
	  printf("%ld %8.4f %9.5f %9.5f %8.3f" ,
		p->index, p->time, p->lon, p->lat, p->depth) ;
	  printf("%2d %2d %2d %6.3f %6.3f %6.3f",
		p->nP, p->nS, p->nIter, p->stdP,p->stdS, p->length ) ; 
	printf("\n") ;
	}
}
int compareDouble( const void *p1, const void *p2 )
{
	double d1,d2 ;
	d1 = *( double *) p1 ;
	d2 = *( double *) p2 ;
	if( d2 < d1 ) return 1 ;
	if( d2 > d1 ) return -1 ;
	return 0 ;
}
double *sortDouble(int n, double *base, int stride, int nPrint )
{
	double *work, *fp ;
	char *cp ;
	int i,n0,n1 ;
	cp = ( char*) base ;
	work = malloc( n * sizeof(double)) ;
	for( i = 0 ; i < n ; i++ ) {
		fp = ( double *) cp ;
		work[i] = *fp ;
		cp += stride ;
	}
	qsort(work,n,sizeof(double),compareDouble) ; 
	if( nPrint == 0 ) return work ;
	n0 = 0 ; n1 = n ;
	if( nPrint > 0 ) { 
		n0 = n - nPrint ;
		if( n0 < 0 ) n0 = 0  ;
	} else {
		n1 = -nPrint ;
		if( n1 > n ) n1 = n ;
	}
	for( i = n0 ; i < n1 ; i++) printf("%8.3f",work[i]) ;
	printf("\n") ;
	return work ;
}
int compareInt( const void *p1, const void *p2 )
{
	return  *(int *) p1 -  *(int *) p2 ;
}
int *sortInt(int n, int *base, int stride, int nPrint )
{
	int *work, *fp ;
	char *cp  ;
	int i,n0,n1 ;
	cp = ( char*) base ;
	work = malloc( n * sizeof(int)) ;
	for( i = 0 ; i < n ; i++ ) {
		fp = ( int *) cp ;
		work[i] = *fp ;
		cp += stride ;
	}
	qsort(work,n,sizeof(int),compareInt) ; 
	n0 = 0 ; n1 = n ;
	if( nPrint > 0 ) { 
		n0 = n - nPrint ;
		if( n0 < 0 ) n0 = 0  ;
	} else {
		n1 = -nPrint ;
		if( n1 > n ) n1 = n ;
	}
	for( i = n0 ; i < n1 ; i++) printf("%6d",work[i]) ;
	printf("\n") ;
	return work ;
}
int compareLL( const void *p1, const void *p2 )
{
	long long *ip1,*ip2 ;
	ip1 = ( long long * ) p1 ;
	ip2 = ( long long * ) p2 ;
	if ( *ip1 > *ip2 ) return 1  ;
	if ( *ip1 < *ip2 ) return -1  ;
	return 0 ;
}
void readIndexList()
{
	int fsize,fd,n,i,j ;
	FILE *ff ;
	fd = open(indexFile,O_RDONLY) ;
	fsize = lseek(fd,0,SEEK_END) ;
	close(fd) ;
	n = 2 + fsize/18 ;
	indexList = malloc(sizeof(long long)*n ) ;
	ff = fopen(indexFile,"r") ;
	for( i = 0 ; i < n ; i++) {
		j = fscanf(ff,"%ld",indexList+i ) ;
/*		printf("%d %d %ld\n",i,j,indexList[i] ) ; */
		if( j < 1 ) break ;
	}
	fclose(ff) ;
	nIndexList = i ;
	printf("nIndexList = %d\n",i) ;
	qsort( indexList,nIndexList, sizeof(long long ) , compareLL ) ;
}
int readTable( char *suffix, int size, void **addr ) 
{
	int fd, space, n, nRec ;
	long long *ip ;
	char fname[100] ;
	sprintf(fname,"%s.%s",baseName,suffix) ;
	fd = open(fname,O_RDONLY )  ;
	space = lseek(fd,0,SEEK_END) ;
	n = lseek(fd,0,SEEK_SET) ;
	*addr = malloc(space+size) ; /* space at end for sentinel */
	n = read(fd,*addr,space) ; 
	nRec = n/size ;
	printf("fname=%s n=%d space=%d size=%d nRec=%d\n",fname,n,space,size,nRec) ;
	if( space != n ) rLog(1,"Error reading %s",fname) ;
	close(fd) ;
	ip = *addr + space ;    *ip = INDEXEND ;
	n = read(fd,*addr,space) ; 
	return nRec ;
}
void printS()
{
	int i ;
	Solution *p ;
	if( shLogLevel < 6 ) return ;
	p = solutions ;
	for( i = 0 ; i < nSol ; i++) {
		printf("%3d %ld %8.3f %8.4f %8.4f %8.1f %3d %3d %3d\n",
			i,p->index,p->time,p->lat,p->lon,p->depth,p->nPhase,p->nP,p->nS  ) ;
		p++ ;
	}
}
void printE()
{
	int i,nn ;
	Event *p ;
	if( shLogLevel < 6 ) return ;
	p = events ;
	nn = nEvent ;
	if( nn > 10 ) nn = 10 ;
	for( i = 0 ; i < nn ; i++) {
		printf("%3d %ld %8.3f %8.4f %8.4f %8.1f\n",
			i,p->index,p->time,p->lat,p->lon,p->depth  ) ;
		p++ ;
	}
}
void printP()
{
	int i,n ;
	Phase *p ;
	if( shLogLevel < 6 ) return ;
	p = phases ;
	n = nPhase ;
	if( n > 500 ) n = 500 ;
	for( i = 0 ; i < n ; i++) {
		printf("%3d %ld %8.3f %s %c\n",i,p->index,p->pTime, p->station,p->type  ) ;
		p++ ;
	}
}
int comparePhase( const void *p1, const void *p2 )
{
	Phase *pp1, *pp2 ;
	long long i1,i2 ;
	pp1 = (Phase *)p1 ; pp2 = (Phase *)p2 ;
	i1 = pp1->index ;
	i2 = pp2->index ;
/*	printf("i1=%ld i2=%ld\n",i1,i2) ; */
	if ( i1 > i2 ) return 1  ;
	if ( i1 < i2 ) return -1  ;
	return pp1->iPhase - pp2->iPhase ;
}
int compareEvent( const void *p1, const void *p2 )
{
	Event *pp1, *pp2 ;
	long long i1,i2 ;
	pp1 = (Event *)p1 ; pp2 = (Event *)p2 ;
	i1 = pp1->index ;
	i2 = pp2->index ;
/*	printf("i1=%ld i2=%ld\n",i1,i2) ; */
	if ( i1 > i2 ) return 1  ;
	if ( i1 < i2 ) return -1  ;
	return 0 ;
}
int testIndexList( long long index )
{
	static int i ;
	int result ;
	if( index == 0 ) { i = 0 ; return 0 ; } /* rewind index list */
	while( indexList[i] < index - 300 ) i++ ;
	if( indexList[i] > index + 600 ) result = 1 ; else result = 0 ;
	if( shLogLevel > 4 ) 
		printf("testIndexList %4d %ld %ld %d %14ld\n",
			i,index,indexList[i],result,index-indexList[i] ) ;
	return result ;
}
int testEvent( Event *ep, Phase *pp, int nP )
{
	int i,ns ;
	static long long lastIndex ;
	long long index ;
	index = ep->index ;
	if(shLogLevel > 4 ) 
		printf(" testEvent : index=%ld  nP=%d\n",index,nP) ;
	if ( index == lastIndex ) return 0 ; /* remove duplicates */
	lastIndex = index ;
	if ( index == (ep+1)->index ) return 0 ;
	if ( ep->lat < latMin ) return 0 ;
	if ( ep->lat > latMax ) return 0 ;
	if ( ep->lon < lonMin ) return 0 ;
	if ( ep->lon > lonMax ) return 0 ;
	if ( ep->depth < zMin ) return 0 ;
	if ( ep->depth > zMax ) return 0 ;
	if( index < indexMin )  return 0 ;
	if( index > indexMax )  return 0 ;
	if( nIndexList) if(testIndexList(index)) return 0 ; 
	if( nP       < phasesMin )  return 0 ; 
	if( nP       > phasesMax )  return 0 ; 
	ns = 0 ; 
	for( i = 0 ; i < nP ; i++) if (pp[i].type == 'S' ) ns++ ;
	if( ns < minSorP ) return 0 ;
	if( nP - ns < minSorP ) return 0 ;
	return 1 ;
}
void countEvents() 
{
	int ie  ;
	long long idx ;
	Phase *pp, *pp1 ;
	Event *ep ;
	pp = phases ;
	ep = events ;
	idx = -1 ;
	for ( ie = 0 ; ie < nEvent ; ie++) {
		if( idx == ep->index )  rLog(1,"duplicate index %ld",(void *) idx) ;
		idx = ep->index ;
		while ( pp->index < idx ) pp++ ;
		pp1 = pp ;
		while ( pp1->index == idx ) pp1++ ;
/*		printf("%5d %ld %3d\n",ie,idx, pp1-pp ) ; */
		nSol += testEvent(ep,pp,pp1-pp) ;
		ep++ ;
		pp = pp1 ;
	}
	printf("Of %d events %d passed\n",nEvent,nSol) ;
}
void checkPhases()
{
	int ip ;
	Phase *pp ;
	for( ip = 0 ; ip < nPhase ; ip++) {
		pp = phases + ip ;
		pp->statP = lookUpStation( pp->station ) ;
	}
}
void checkPhasesX() /* OBSOLETE remove phases to far from source. Distance estimated from phase times.*/
{	
	int ie,ip ;
	double ttime,dist ;
	long long idx  ;
	Phase *pp,*pp1 ;
	Event *ep ;
	pp = phases ;
	ep = events ;
	for( ie = 0 ; ie < nEvent ; ie++ ) {
		idx = ep->index ;
		while ( pp->index < idx ) pp++ ;
		pp1 = pp ; 
		while ( pp1->index == idx ){
			ttime = pp1->pTime - ep->time ;
			if ( pp1->type == 'S' ) dist = ttime*3.6 ;
			if ( pp1->type == 'P' ) dist = ttime*6.4 ;
			if( dist > distMax ) pp1->type = 'x' ;
/*			printf("%5d %ld %ld %8.3f %8.3f %8.3f %8.3f %c\n",
				ie,ep->index,pp1->index,ep->time,pp1->pTime,ttime,dist,pp1->type) ; */
			pp1++ ;
		}
		ep++ ;
		pp = pp1 ;
	}
	pp = phases ;
	pp1 = phases ;
	for(ip = 0 ; ip < nPhase ; ip++) {
		*pp1 = *pp ;
		if (pp->type != 'x' ) pp1++ ;
		*pp++ ;
	}
	printf("%d phases further than %8.2f km from source\n",pp -pp1,distMax) ;
	nPhase = pp1 - phases ;
	for( ip = 0 ; ip < nPhase ; ip++) {
		pp = phases + ip ;
		pp->statP = lookUpStation( pp->station ) ;
	}
}
void makeSolutions()
{
	int ie  ;
	long long idx ;
	Phase *pp, *pp1 ;
	Event *ep ;
	Solution *sp ;
	pp = phases ;
	ep = events ;
	solutions = calloc(nSol,sizeof(Solution)) ;
	sp = solutions ;
	for ( ie = 0 ; ie < nEvent ; ie++) {
		idx = ep->index ;
		while ( pp->index < idx ) pp++ ;
		pp1 = pp ;
		while ( pp1->index == idx ) pp1++ ;
		if( testEvent(ep,pp,pp1-pp)) {
			sp->index = ep->index ;
			sp->lat   = ep->lat ;
			sp->lon   = ep->lon ;
			sp->depth = ep->depth ;
			sp->time  = ep->time ;
			sp->phase = pp ;
			sp->nPhase = pp1 - pp ;
			sp++ ;
		}
		ep++ ;
		pp = pp1 ;
	}
/*	free(events) ; */
}
initVel()
{
	char pv[60],sv[60] ;
	sprintf(pv,"%sp.vel",velPrefix ) ;
	sprintf(sv,"%ss.vel",velPrefix ) ;
	initVelModel(20,&mp) ; readVelModel(pv,&mp) ; 
	initVelModel(20,&ms) ; readVelModel(sv,&ms) ; 
	initVelModel(20,&initP) ; readVelModel(pv,&initP) ; 
	initVelModel(20,&initS) ; readVelModel(sv,&initS) ; 
	
/*	vFInitFromMemory() ; */
}
void pass1()
{	Solution *lp, *op ;
	Phase *pp ;
	int ii,sumi,nn ;
/*	vFInitFromMemory() ; */
	sumi = 0 ;
	printf("enter pass1: nSol = %d\n",nSol);
	for( ii = 0 ; ii < nSol ; ii++) {
		lp = solutions + ii ;
		lp->nIter = 1 ;
		pp = lp->phase ;
		locate( lp,pp) ;	
		if(shLogLevel > 0 )
		   printf("time=%8.3f ",lp->time ) ;
		   printf("nP=%3d nS=%3d nIter=%3d stdP=%10.4f stdS=%10.4f length=%10.4f\n",
			lp->nP,lp->nS,lp->nIter,lp->stdP, lp->stdS,lp->length) ;
		sumi += lp->nIter ;
	}
	op = solutions ;
	lp = solutions ;
	for( ii = 0 ; ii < nSol ; ii++) {
		if( 	(lp->nIter < 12 ) &&
			(lp->stdP < 0.07 ) &&
			(lp->stdS < 0.1 )) *op++ = *lp ;
		lp++ ;
	}
	printf("%d iterations, nSol = %d nWorking = %d\n",sumi, nSol,op-solutions) ; 
	nSol = op - solutions ; 
	printf("leave pass1: nSol = %d\n",nSol);

}
double processVel() 
{
	Solution work, *lp, *op ;
	Phase  *pp ;
	double sumP,sumS,yy ;
	int nP,nS ;
	int i,n ;
	nP = 0 ; nS = 0 ;
	sumP = 0.0 ; sumS = 0.0 ;
	n = nSol ;
	op = solutions ;
	for ( i = 0 ; i < n ; i++ ) {
		lp = solutions + i ;
		work = *lp ;
		pp = work.phase ;
		locate( &work,pp) ;
		nP += work.nP ;
		nS += work.nS ;
		sumP += work.sumP ;
		sumS += work.sumS ;
		*op++ = work ;
	}
	nSol = op-solutions ;
	yy = ( sumP*(1.0-weightS) + sumS*weightS ) / ( nP*(1.0-weightS) + nS*weightS ) ;
	yy = sqrt(yy) ;
	return 1000*yy ; /* convert seconds to milliseconds */
}
double processVelOld( Solution *sol, int n )
{
	Solution *lp ;
	Phase *pp ;
	int i ;
	double sst ;
	sst = 0 ;
/*	vFInitFromMemory() ;  */
	for( i = 0 ; i < n ; i++ ) {
		lp = sol + i ;
		pp = lp->phase ;
		locate( lp,pp) ;
		sst +=  ( 0.5 * lp->stdP + 0.3 * lp->stdS ) ;
	}
	(void) sortInt(nSol,&(solutions->nIter),sizeof(Solution),10 );
	return sst/n ;
}
double lowerParLimit( ipar )
{		/* find lowest value of velocity that avoids increase in gradient */
	int i ;
	double slope, limit ;
	VelModel *m ;
	if ( ipar % 2 ) m = &ms ;
		else m = &mp ;
	i = skipLayers + ipar/2 ;
	if( i == 0 ) return 0.0 ;
	slope = ( m->v[i+1] - m->v[i-1] ) / ( m->z[i+1] - m->z[i-1] ) ;
	limit = m->v[i-1] + ( m->z[i] - m->z[i-1] ) * slope ;
	return limit ;
}
double upperParLimit( ipar )
{		/* find higest value of velocity that avoids increase in gradient */
	int i ;
	double slope, limit1, limit2 ;
	VelModel *m ;
	if ( ipar % 2 ) m = &ms ;
		else m = &mp ;
	i = skipLayers + ipar/2 ;
	slope = ( m->v[i+2] - m->v[i+1] ) / ( m->z[i+2] - m->z[i+1] ) ;
	limit2 = m->v[i+1] - ( m->z[i+1] - m->z[i] ) * slope ;
	if( i < 2 ) return limit2 ;
	slope = ( m->v[i-1] - m->v[i-2] ) / ( m->z[i-1] - m->z[i-2] ) ;
	limit1 = m->v[i-1] + ( m->z[i] - m->z[i-1] ) * slope ;
	if( limit1 > limit2 ) return limit2 ;
	return limit1 ;
}
void testLimits(int nvel)
{
	int i ;
	double p1,p2,s1,s2,slopeP,slopeS ;
	for( i = 0 ; i < nvel-skipLayers ; i++ ) {
		p1 = lowerParLimit(2*i) ;
		p2 = upperParLimit(2*i) ;
		s1 = lowerParLimit(2*i+1) ;
		s2 = upperParLimit(2*i+1) ;
		slopeP = (mp.v[i+1] - mp.v[i]) / ( mp.z[i+1]-mp.z[i] ) ;
		slopeS = (ms.v[i+1] - ms.v[i]) / ( ms.z[i+1]-ms.z[i] ) ;
		printf("%8.2f %8.2f %8.2f %8.2f %8.3f %8.2f %8.2f %8.2f %8.3f\n",
			mp.z[i],p1,mp.v[i],p2,slopeP,s1,ms.v[i],s2,slopeS) ;
	}
}
int metropolisTest( double y1, double y0 ) 
{
	double t, rr ;
	if( y1 < y0 ) return 1 ;
	if( 0 == metropolis ) return 0 ;
	t = exp( - ( y1 - y0 ) / temperature ) ;
	rr = random() * 1.0 / RAND_MAX ;
	printf("y1=%9.6f y0=%9.6f t=%9.5f rr=%9.5f temperature=%10.8f\n",y1,y0,t,rr,temperature) ;
	temperature *= 0.999 ;
	if (t > rr) return 1 ;
	return 0 ;
}
void searchRandom( int nVel )
{
	int i, nPar, iPar, nOk, lastNSol ;
	double *value[50],**vp, work,y0,y1  ;
	double range,delta,rangeScale ;
	double upperL, lowerL ;
	time_t tt ;
	
	vp = value ;
	for( i = skipLayers ; i < nVel ; i++){
		*vp++ = mp.v+i ;
		*vp++ = ms.v+i ;
	}
	nVel = nLayers - skipLayers ;
	i = 0 ;
	if( noRandomize == 0 ) srandom(time(&tt)) ;
	nOk = 0 ;
	nPar = 2*nVel ;
	range = 0.060 ;
	y0 = processVel() ; 
	do {
		iPar = shuffle(nPar) ;
		iPar %= nPar ;
		upperL = upperParLimit(iPar) ;
		lowerL = lowerParLimit(iPar) ;
		work = *value[iPar] ; 
		if( gaussFlag) {
			*value[iPar] = grandom(lowerL,upperL,work,range) ;
		} else {
			if( (work + range) < upperL ) upperL = work + range ;
			if( (work - range) > lowerL ) lowerL = work - range ;
			rangeScale = RAND_MAX/(upperL - lowerL) ;
			*value[iPar] = lowerL + (random() / rangeScale) ;
		}
		delta = *value[iPar] - work ;
/*		printf("i=%d iPar= %d delta=%8.2f \n",i,iPar,delta) ; */
		lastNSol = nSol ;
		y1 = processVel() ; 
		if(( metropolisTest( y1 , y0 ) ) && ( nSol == lastNSol )) {  
			nOk++ ;
			y0 =  y1 ;
		} else {
			*value[iPar] = work ;
			if( nSol != lastNSol ) 
				y0 = processVel() ;
			range *= 0.999 ;
		}
		printf("%3d ipar=%2d range=%9.5f delta=%9.5f y0=%10.7f %3d %3d %4d\n",
			i,iPar,range,delta,y0,nOk,i-nOk,nSol ) ;
	} while (i++ < nIterations ) ;
	for( i = 0 ; i < nPar ; i++) printf("%8.4f", *value[i] ) ;
	printf("\n") ;
}
void search( int nVel )
{
	int i,nPar,iPass ;
	double *value[50], **vp, *work,x1,x2,x3,y1,y2,y3  ;	
	double xstep, damper,limit ;
	double a,b,c,xx ;
	damper = 0.3 ;
	limit = 0.2 ;
	xstep = 0.04 ;
	vp = value ;
	for (i = 0; i < nVel ; i++ ) {
		*vp++ = mp.v+i	;
		*vp++ = ms.v+i	;
	}
	nPar = 2*nVel ;
	for( iPass = 0 ; iPass < 10 ; iPass++ ) {
	    for( i = 0 ; i < nPar ; i++) {
		work = value[i] ;
		x2 = *work ;
		y2 = processVel(solutions,nSol) ;
		x1 = x2 - xstep ;
		*work = x1 ;
		y1 = processVel(solutions,nSol) ;
		x3 = x2 + xstep ;
		*work = x3 ;
		y3 = processVel(solutions,nSol) ;
		c = 0.5*(y3+y1) - y2 ;
		b = 0.5*(y3-y1) ;
		a = y2  ;
		xx = -0.5*b*xstep/c ;
		if( xx > limit ) xx = limit ;
		if( xx < -limit ) xx = -limit ;
		*work = x2 + damper * xx ;
	        printf("%9.5f%9.5f%9.5f |",y2,xx,c) ;
	     }
	     printf("\n") ;
	}
	printVelModel(&mp) ;
}
void getData()
{
	nEvent = readTable("event", sizeof(Event),(void *) &events) ;
	nPhase = readTable("phase", sizeof(Phase),(void *) &phases) ;
	printE() ;
/*	printP() ;
	printE() ; */
	qsort(phases,nPhase,sizeof(Phase), comparePhase ) ; 
	qsort(events,nEvent,sizeof(Event), compareEvent ) ; 
	printP() ;
/*	printP() ;
	printE() ; */
	printf("%d events, %d phases\n",nEvent,nPhase) ;
	if( indexFile ) readIndexList() ;
	checkPhases() ;
	countEvents() ;
	testIndexList(0) ; /* rewind index list */
	printf("nEvent=%d nSol=%d\n",nEvent,nSol) ;
	printE() ;
	makeSolutions() ;
	printS() ;
	printP() ;
}
void doIt()
{
	getData() ;
/*	printP() ; */
	initVel() ;
	pass1() ;
/*	serach(6) ;    */
	searchRandom(nLayers) ;
}
void doLocations()
{
	int i ;
	Solution *lp ;
	Phase *pp ;
	getData() ;
	initVel() ;
	pass1() ;
/*	if( rayTrace == 0 )  vFInitFromMemory() ;  */
	for(i = 0  ; i < nSol ; i++ ) {
		lp = solutions + i ;
		pp = lp->phase ;
		locate(lp,pp) ;
	}
}
void printReport(int ac, char **av)
{
	int j ;
	time_t tt ;
/*	(void) sortDouble(nSol,&(solutions->lat),sizeof(Solution) ); */
	printf("nIter : ") ;
	(void) sortInt(nSol,&(solutions->nIter),sizeof(Solution),10 );
	printf("nP    : ") ;
	(void) sortInt(nSol,&(solutions->nP),sizeof(Solution),10 );
	printf("nP -  : ") ;
	(void) sortInt(nSol,&(solutions->nP),sizeof(Solution),-10 );
	printf("nS    : ") ;
	(void) sortInt(nSol,&(solutions->nS),sizeof(Solution),10 );
	printf("nS -  : ") ;
	(void) sortInt(nSol,&(solutions->nS),sizeof(Solution),-10 );
	printf("length: ") ;
	(void) sortDouble(nSol,&(solutions->length),sizeof(Solution), 10) ;
	printf("depth : ") ;
	(void) sortDouble(nSol,&(solutions->depth),sizeof(Solution),  5 ) ;
	printf("depth-: ") ;
	(void) sortDouble(nSol,&(solutions->depth),sizeof(Solution), -5 ) ;
	printParameters(stdout) ;
	for(j = 0 ; j < ac ; j++ ) printf("%s ",av[j] ) ; 
	tt = time(&tt) ;
	printf("\n%s",ctime(&tt)) ;
}
void printPhaseTable()
{
	Phase *pp ;
	Station *sp ;
	int i ;
	nPhase = readTable("phase", sizeof(Phase),(void *) &phases) ;
	printf("\nIn printPhaseTable: nPhase = %d\n",nPhase) ;
	qsort(phases,nPhase,sizeof(Phase), comparePhase ) ;
	for( i = 0 ; i < nPhase ; i++ ) {
		pp = phases + i ;
		 sp = pp->statP ; 
		printf("%ld %s %c %8.3f\n",
			pp->index,pp->station,pp->type, pp->pTime)  ; 
	}
}
void printIndexTable()
{
	Event *ep ;
	static struct tm tm ;
	int i ;
	long long tw ;
	time_t ttIndex ;
	char *cp ;
	nEvent = readTable("event", sizeof(Event),(void *) &events) ;
	qsort(events,nEvent,sizeof(Event), compareEvent ) ; 
	printf("printIndexTable: nEvent = %d\n",nEvent) ;
	for( i = 0 ; i < nEvent ; i++) {
		ep = events+i ;
		tw = ep->index / 1000 ;
		tm.tm_sec  = tw%100 ;	tw /= 100 ;
		tm.tm_min  = tw%100 ;	tw /= 100 ;
		tm.tm_hour = tw%100 ;	tw /= 100 ;
		tm.tm_mday  = tw%100 ;	tw /= 100 ;
		tm.tm_mon  = tw%100 ;	tw /= 100 ;
		tm.tm_year = tw - 1900 ;	
		ttIndex = mktime(&tm) ;
		ttIndex += floor(ep->time) ;
		tm = *gmtime(&ttIndex) ;
		cp = asctime(&tm) ;
		printf("%ld %ld %9.3f %s",ep->index,ep->index/1000,ep->time,cp) ;
	}
	printPhaseTable() ;
	exit(0) ;
}
void printResults(char *resultFile) 
{
	FILE *ff ;
	Solution *p ;
	int i ;
	int nP,nS ;
	double sumP,sumS ;
	ff = fopen(resultFile,"w") ;
	for( i = 0 ; i < nLayers ; i++) {
		fprintf(ff,"RV %8.3f %8.3f %8.3f %8.3f %8.3f\n",
			initP.z[i],initP.v[i],mp.v[i],initS.v[i],ms.v[i] ) ;
	}
	sumP = 0.0 ; sumS = 0.0 ; nP = 0 ; nS = 0 ;
	for( i = 0 ; i < nSol ; i++) {
	  p = solutions + i ;
	  nP += p->nP ;
	  nS += p->nS ;
	  sumP += p->sumP ;
	  sumS += p->sumS ;
	  fprintf(ff,"RS %ld %8.4f %9.5f %9.5f %8.3f" ,
		p->index, p->time, p->lon, p->lat, p->depth) ;
	  fprintf(ff,"%3d %3d %3d %6.3f %6.3f %6.3f",
		p->nP, p->nS, p->nIter, p->stdP,p->stdS, p->length ) ; 
	  fprintf(ff,"\n") ;
	}
	fclose(ff) ;
	
}
void testing()
{
	testLimits(8) ;
}
int main( int ac, char **av)
{
int cc ;
	extern char *optarg ;	
	feenableexcept(FE_INVALID) ;
	shLogLevel = 2 ;
	rayTrace = 1 ;
	while( EOF != ( cc = getopt(ac,av,"aAgtwd:e:z:Z:b:B:l:L:n:N:m:i:I:r:R:D:W:vspT:M:j:V:XS:P:"))) {
		switch(cc) {
		case 'a' : metropolis = 1 ; break ;
		case 'A' : noRandomize = 1 ; break ;
		case 'g' : gaussFlag = 1 ; break ;  /* gaussian distribution in model space */
		case 'S' : nIterations = atoi(optarg) ; break ;
		case 'd' : shLogLevel = atoi(optarg) ; break ;
		case 'w' : locateSILWeight = 1 ; break ;
/*		case 't' : rayTrace = 0 ; break ; */
		case 'e' : baseName = optarg ; break ;
		case 'z' : zMin = atof(optarg) ;  break ;
		case 'Z' : zMax = atof(optarg) ;  break ;
		case 'b' : latMin = atof(optarg) ; break ;
		case 'B' : latMax = atof(optarg) ; break ;
		case 'l' : lonMin = atof(optarg) ; break ;
		case 'L' : lonMax = atof(optarg) ; break ;
		case 'n' : phasesMin = atoi(optarg) ; break ;
		case 'N' : phasesMax = atoi(optarg) ; break ;
		case 'm' : minSorP = atoi(optarg) ; break ;
		case 'i' : indexMin = atoll(optarg) ; 
			/* appending zeroes to index limits as needed */
			while (indexMin < INDEXEND/10 ) indexMin *= 10 ; 
			break ;
		case 'I' : indexMax = atoll(optarg)  ;
			while (indexMax < INDEXEND/10 ) indexMax *= 10 ;
			; break ;
		case 'r' : indexFile = optarg ; break ;
		case 'R' : indexFileFraction = atoi(optarg) ; break ;
		case 'D' : distMax = atof(optarg) ; break ;
		case 'W' : weightS = atof(optarg) ; break ;
		case 'T' : temperature = atof(optarg) ; break ;
		case 'M' : nLayers = atoi(optarg) ; break ;
		case 'j' : skipLayers = atoi(optarg) ; break ;
		case 'V' : velPrefix = optarg ; break ;
/* action flags  */
		case 'v' : doIt() ; break ;
		case 's' : doLocations() ; break ;
		case 'p' : printSolutions() ; break ;
		case 'X' : printIndexTable();  break ;
		case 'P' : printResults(optarg); break ;
	}}
	printReport(ac,av) ;
	return 0 ;
}

