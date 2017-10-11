
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

char *baseName = "testing" ;
Event *events ;
Phase *phases ;
Solution *solutions ;
int nEvent, nPhase, nSol ;
extern VelModel mp,ms ;

/*		Default values */
int phasesMin = 8 ;
int phasesMax = 999 ;
int minSorP = 2 ;
double zMax = 8.0 ;
double zMin = 2.5 ; 
double latMin = 63.0 ;
double latMax = 65.0 ;
double lonMin = -25.0 ;
double lonMax = -18.0 ; 
double distMax  = 90.0 ;
long long indexMin = 199107 ;
long long indexMax = INDEXEND ;

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
void printE()
{
	int i ;
	Event *p ;
	p = events ;
	for( i = 0 ; i < 8 ; i++) {
		printf(" %ld %8.3f\n",p->index,p->time  ) ;
		p++ ;
	}
}
void printP()
{
	int i,n ;
	Phase *p ;
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
int testEvent( Event *ep, Phase *pp, int nP )
{
	int i,ns ;
	long long index ;
	index = ep->index ;
	if ( index == (ep-1)->index ) return 0 ; /* remove duplicates */
	if ( index == (ep+1)->index ) return 0 ;
	if ( ep->lat < latMin ) return 0 ;
	if ( ep->lat > latMax ) return 0 ;
	if ( ep->lon < lonMin ) return 0 ;
	if ( ep->lon > lonMax ) return 0 ;
	if ( ep->depth < zMin ) return 0 ;
	if ( ep->depth > zMax ) return 0 ;
	if( index < indexMin )  return 0 ;
	if( index > indexMax )  return 0 ;
	if( nP       < phasesMin )  return 0 ; 
	if( nP       > phasesMax)  return 0 ; 
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
void checkPhases() /* remove phases to far from source. Distance estimated from phase times.*/
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
	printf("%d phases further than %f km from source\n",pp -pp1,distMax) ;
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
	free(events) ;
}
initVel()
{
	initVelModel(20,&mp) ; readVelModel("silp.vel",&mp) ; 
	initVelModel(20,&ms) ; readVelModel("sils.vel",&ms) ; 
	
/*	vFInitFromMemory() ; */
}
void pass1()
{	Solution *lp ;
	Phase *pp ;
	int ii,sumi ;
	for( ii = 0 ; ii < nSol ; ii++) {
		lp = solutions + ii ;
		lp->nIter = 1 ;
		pp = lp->phase ;
		locate( lp,pp) ;	
		if(shLogLevel > 3 )
		   printf("nP=%3d nS=%3d nIter=%3d stdP=%10.4f stdS=%10.4f length=%10.4f\n",
			lp->nP,lp->nS,lp->nIter,lp->stdP, lp->stdS,lp->length) ;
		sumi += lp->nIter ;
	}

}
void doIt()
{
	nEvent = readTable("event", sizeof(Event),(void *) &events) ;
	events[nEvent].index = INDEXEND ;
	nPhase = readTable("phase", sizeof(Phase),(void *) &phases) ;
	phases[nPhase].index = INDEXEND ;
/*	printP() ;
	printE() ; */
	qsort(phases,nPhase,sizeof(Phase), comparePhase ) ; 
	qsort(events,nEvent,sizeof(Event), comparePhase ) ; 
/*	printP() ;
	printE() ; */
	printf("%d events, %d phases\n",nEvent,nPhase) ;
	checkPhases() ;
	countEvents() ;
	makeSolutions() ;
/*	printP() ; */
	initVel() ;
	pass1() ;
}

int main( int ac, char **av)
{
int cc ;
	extern char *optarg ;	
	feenableexcept(FE_INVALID) ;
	shLogLevel = 4 ;
	while( EOF != ( cc = getopt(ac,av,"e:z:Z:b:B:l:L:n:N:m:i:I:D:"))) {
		switch(cc) {
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
		case 'i' : indexMin = atoll(optarg) ; break ;
		case 'I' : indexMax = atoll(optarg) ; break ;
		case 'D' : distMax = atof(optarg) ; break ;
	}}
	while (indexMin < INDEXEND/10 ) indexMin *= 10 ; /* appending zeroes to index limits as neede */
	while (indexMax < INDEXEND/10 ) indexMax *= 10 ;
	printf("Index range: %ld %ld Phase range: %d %d\n", indexMin, indexMax, phasesMin,phasesMax ) ;
	doIt() ;
	return 0 ;
}

