/*  read sil data from .eve files */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "../rayt/ray.h"

char *outputBaseName = "testing" ;
char *rootName = "/data/bergth/4D/testing/sk1/1998/apr" ;

int month2number( char *monthName)
{
	static char *mn[12] = { "jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec" } ;
	int i ;
	for( i = 0 ; i < 12 ; i++) {
		if( 0==strcmp( monthName, mn[i] ) ) return(i) ;
	}
	return -1 ;
}
writePhase( Phase *pp )
{
	static FILE *ff ;
	char fname[200] ;
	if( pp == NULL ) {
		fclose(ff) ;
		return ;
	}
	if( NULL == ff ) {
		sprintf(fname,"%s.phase",outputBaseName) ;
		ff = fopen(fname,"w") ;
	}
	fwrite(pp, sizeof(Phase),1,ff) ;
}
writeEvent( Event *pp )
{
	static FILE *ff ;
	char fname[200] ;
	if( pp == NULL ) {
		fclose(ff) ;
		return ;
	}
	if( NULL == ff ) {
		sprintf(fname,"%s.event",outputBaseName) ;
		ff = fopen(fname,"w") ;
	}
	fwrite(pp, sizeof(Event),1,ff) ;
}
void readEFile( char *ename ) 
{
	FILE *efile ;
	char *bp, *cp ;
	static struct tm tm ;
	time_t ttIndex, ttEvent,ttPhase ;
	int indexms;
	static size_t n=100 ;
	int len,j,iPhase ;
	char line[100], buffer[100] ;
	double eventSec, sourceTime,phaseSec,eventFrac,phaseFrac,dist,ttime ;
	static Event ev ;
	static Phase phase ;
	Station *sp ;
	cp=strncpy(buffer,ename,100) ;
	if( shLogLevel > 4 )  printf("%s",ename) ; 
	while(*cp){ 
		if ((*cp == '/') ||(*cp == ':')) {*cp = 0 ; bp = cp+1 ; }
		cp++ ;
	}
	ename[cp-buffer-1] = 0 ;
	indexms = atoi(bp) ;
	efile = fopen(ename,"r") ;
	cp = fgets(line,n,efile) ;
	j = atoi(line+60) ;
	tm.tm_sec = j%100 ;
	j = j / 100 ;
	tm.tm_min = j%100 ;
	tm.tm_hour = j/100 ;
	tm.tm_mday=atoi(bp-18) ;
	tm.tm_mon=month2number(bp-22) ;
	tm.tm_year=atoi(bp-27)-1900 ;
	ev.index = (((1900+tm.tm_year)*100+1+tm.tm_mon)*100+tm.tm_mday)*100+tm.tm_hour ;
	ev.index = ((ev.index*100+tm.tm_min)*100+tm.tm_sec)*1000+indexms ; 
	ttIndex = mktime(&tm) ;
	cp = asctime(&tm) ;
	if( shLogLevel > 4 ) {
		printf("%ld ",ev.index) ;
		printf("%s %s ",cp,bp);
		printf(" ename =%s_\n",ename) ;
	}
	cp = fgets(line,n,efile) ;
	tm.tm_year = atoi(line+12) ;
	tm.tm_mon = atoi(line+15) - 1 ;
	tm.tm_mday= atoi(line+18) ;
	tm.tm_hour= atoi(line+23) ;
	tm.tm_min = atoi(line+27) ;
	eventSec = atof(line+31) ;
	tm.tm_sec = trunc(eventSec) ; 
	eventFrac = eventSec - tm.tm_sec ;
	ttEvent = mktime(&tm) ;
	sourceTime = (ttEvent-ttIndex) + eventSec ;
	if( shLogLevel > 4 ) 
		printf("eventSec = %f sourceTime = %f \n",eventSec,sourceTime) ;
	cp = fgets(line,n,efile) ; ev.lat = atof(line+10) ;
	cp = fgets(line,n,efile) ; ev.lon = atof(line+10) ;
	fgets(line,n,efile) ; ev.depth = atof(line+12) ;
	ev.time = eventFrac + ttEvent - ttIndex ;
	if( shLogLevel > 4 ) 
		printf("lat lon depth time: %f %f %f %f\n",ev.lat,ev.lon,ev.depth,ev.time) ;
 	writeEvent(&ev) ; 
	cp = fgets(line,n,efile) ;
	cp = fgets(line,n,efile) ;
	iPhase = 0 ;
	while (*cp == ' ' ) {
		tm.tm_hour = atoi(line+6) ;
		tm.tm_min  = atoi(line+9) ;
		phaseSec = atof(line+12);
		tm.tm_sec =  trunc(phaseSec) ;
		phaseFrac = phaseSec - tm.tm_sec ;
		ttPhase = mktime(&tm) ;
		dist = atof(line+31) ;
		ttime = ttPhase-ttEvent+phaseFrac-eventFrac ;
		line[4] = 0 ;
		sp = lookUpStation(line+1) ;
		phase.index = ev.index ;
		phase.pTime = ttPhase-ttIndex + phaseFrac ;
		phase.type = line[5] ;
		phase.iPhase = iPhase++ ;
		strncpy(phase.station,line+1,3) ;
		phase.station[3] = 0 ;
		if( shLogLevel > 4 ) 
			printf(" hour,min sec ttime %d %d %f %d %f %f %f %s\n",
			tm.tm_hour,tm.tm_min,phaseSec,ttPhase-ttEvent,dist,ttime,dist/ttime,sp->name) ;
		writePhase(&phase) ;
		cp = fgets(line,n,efile) ;
	}
	fclose(efile) ;
}
void getEList(  )
{
	FILE *efiles ;
	char cmd[100] ;
	char *line ;	
	size_t n ;
	int len ;
	sprintf(cmd,"find %s -name *.eve ",rootName) ;
	printf("%s\n",cmd) ;
	efiles = popen(cmd,"r") ;
	while ( 0 <   (len = getline(&line,&n,efiles))) { 
		readEFile(line) ;
/*		printf("%d %d _ %s",len, n,line) ; */
	} while (len > 0) ;
	writeEvent(NULL) ; /* close output files */
	writePhase(NULL) ;
	
}
int main(int ac , char **av)
{
	int cc ;
	extern char *optarg ;
	shLogLevel = 2 ;
	while( EOF != ( cc = getopt(ac,av,"l:b:o:"))) {
		switch(cc) {
		case 'l' : shLogLevel = atoi(optarg) ; break ;
		case 'b' : rootName = optarg ; break ;
		case 'o' : outputBaseName = optarg ; break ;
	}}
	getEList();
}
