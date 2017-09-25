
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
int nEvent, nPhase ;

int readTable( char *suffix, int size, void **addr ) 
{
	int fd, space, n, nRec ;
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
doIt()
{
	nEvent = readTable("event", sizeof(Event),(void *) &events) ;
	events[nEvent].index = INDEXEND ;
	nPhase = readTable("phase", sizeof(Phase),(void *) &phases) ;
	phases[nPhase].index = INDEXEND ;
	printP() ;
	printE() ;
	qsort(phases,nPhase,sizeof(Phase), comparePhase ) ; 
	qsort(events,nEvent,sizeof(Event), comparePhase ) ; 
	printP() ;
	printE() ;
	
}

int main( int ac, char **av)
{
	doIt() ;
}

