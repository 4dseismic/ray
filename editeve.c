
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "../rayt/ray.h"

char *baseName = "testing" ;
Event *events ;
Phase *phases ;
Solution *solutions ;
int nEvent, nPhase, nSol ;
/*		Default values */
int phasesMin = 5 ;
int phasesMax = 999 ;
double zMax = 8.0 ;
double zMin = 2.5 ; 
double latMin = 63.0 ;
double latMax = 65.0 ;
double lonMin = -25.0 ;
double lonMax = -18.0 ; 
double distMax  = 90.0 ;
long long indexMin = 1900 ;
long long indexMax = INDEXEND;

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
	int i ;
	Phase *p ;
	p = phases ;
	for( i = 0 ; i < 18 ; i++) {
		printf(" %ld %8.3f %s %c\n",p->index,p->pTime, p->station,p->type  ) ;
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
	if ( ep->lat < latMin ) return 0 ;
	if ( ep->lat > latMax ) return 0 ;
	if ( ep->lon < lonMin ) return 0 ;
	if ( ep->lon > lonMax ) return 0 ;
	if ( ep->depth < zMin ) return 0 ;
	if ( ep->depth > zMax ) return 0 ;
	if( ep->index < indexMin )  return 0 ;
	if( ep->index > indexMax )  return 0 ;
	if( nP       < phasesMin )  return 0 ; 
	if( nP       > phasesMax)  return 0 ; 
	return 1 ;
}
void checkEvents() 
{
	int ie  ;
	long long idx ;
	Phase *pp, *pp1 ;
	Event *ep ;
	pp = phases ;
	ep = events ;
	for ( ie = 0 ; ie < nEvent ; ie++) {
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
	double ts,tp,ttime ;
	long long idx  ;
	Phase *pp,*pp1 ;
	Event *ep ;
	ts = distMax/3.5 ;
	tp = distMax/6.4 ;
	pp = phases ;
	ep = events ;
	for( ie = 0 ; ie < nEvent ; ie++ ) {
		idx = ep->index ;
		while ( pp->index < idx ) pp++ ;
		pp1 = pp ; 
		while ( pp1->index == idx ){
			ttime = pp1->pTime - ep->time ;
			if ( pp1->type == 's' ) 
				if( ttime > ts ) pp1->type = 'x' ;
			if ( pp1->type == 'p' ) 
				if( ttime > tp ) pp1->type = 'x' ;
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
	checkEvents() ;
	makeSolutions() ;
}

int main( int ac, char **av)
{
int cc ;
	extern char *optarg ;	
	while( EOF != ( cc = getopt(ac,av,"e:z:Z:b:B:l:L:n:N:i:I:D:"))) {
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
		case 'i' : indexMin = atoll(optarg) ; break ;
		case 'I' : indexMax = atoll(optarg) ; break ;
		case 'D' : distMax = atof(optarg) ; break ;
	}}
	while (indexMin < INDEXEND/10 ) indexMin *= 10 ;
	while (indexMax < INDEXEND/10 ) indexMax *= 10 ;
	printf("Index range: %ld %ld Phase range: %d %d\n", indexMin, indexMax, phasesMin,phasesMax ) ;
	doIt() ;
	return 1 ;
}

