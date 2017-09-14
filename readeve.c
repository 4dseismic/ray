/*  read sil data from .eve files */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "../rayt/ray.h"

int month2number( char *monthName)
{
	static char *mn[12] = { "jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec" } ;
	int i ;
	for( i = 0 ; i < 12 ; i++) {
		if( 0==strcmp( monthName, mn[i] ) ) return(i) ;
	}
	return -1 ;
}
void readEFile( char *ename ) 
{
	FILE *efile ;
	char *bp, *cp ;
	static struct tm tm ;
	time_t tt, ttEvent ;
	int indexms;
	static size_t n=100 ;
	int len ;
	char line[100], buffer[100] ;
	double fsec, sourceTime ;
	static Solution sol ;
	cp=strncpy(buffer,ename,100) ;
	printf("%s",ename) ; 
	while(*cp){ 
		if ((*cp == '/') ||(*cp == ':')) {*cp = 0 ; bp = cp+1 ; }
		cp++ ;
	}
	ename[cp-buffer-1] = 0 ;
	indexms = atoi(bp) ;
	tm.tm_sec=atoi(bp-3) ;
	tm.tm_min=atoi(bp-6) ;
	tm.tm_hour=atoi(bp-15) ;
	tm.tm_mday=atoi(bp-18) ;
	tm.tm_mon=month2number(bp-22) ;
	tm.tm_year=atoi(bp-27)-1900 ;
	sol.index = (((1900+tm.tm_year)*100+1+tm.tm_mon)*100+tm.tm_mday)*100+tm.tm_hour ;
	sol.index = ((sol.index*100+tm.tm_min)*100+tm.tm_sec)*1000+indexms ; 
	tt = mktime(&tm) ;
	cp = asctime(&tm) ;
	printf("%ld ",sol.index) ;
	printf("%s %s ",cp,bp);
	printf(" ename =%s_\n",ename) ;
	efile = fopen(ename,"r") ;
	cp = fgets(line,n,efile) ;
	cp = fgets(line,n,efile) ;
	tm.tm_year = atoi(line+12) ;
	tm.tm_mon = atoi(line+15) - 1 ;
	tm.tm_mday= atoi(line+18) ;
	tm.tm_hour= atoi(line+23) ;
	tm.tm_min = atoi(line+27) ;
	tm.tm_sec = 0 ;
	ttEvent = mktime(&tm) ;
	fsec = atof(line+31) ;
	sourceTime = (ttEvent-tt) + fsec ;
	printf("fsec = %f sourceTime = %f \n",fsec,sourceTime) ;
	cp = fgets(line,n,efile) ; sol.lat = atof(line+10) ;
	cp = fgets(line,n,efile) ; sol.lon = atof(line+10) ;
	fgets(line,n,efile) ; sol.depth = atof(line+12) ;
	printf("lat lon depth: %f %f %f\n",sol.lat,sol.lon,sol.depth) ;
	cp = fgets(line,n,efile) ;
	cp = fgets(line,n,efile) ;
	while (*cp != '0' ) {
		printf(cp) ;
		cp = fgets(line,n,efile) ;
	}
}
void getEList( char *rootname )
{
	FILE *efiles ;
	char cmd[100] ;
	char *line ;	
	size_t n ;
	int len ;
	sprintf(cmd,"find %s -name *.eve ",rootname) ;
	printf("%s\n",cmd) ;
	efiles = popen(cmd,"r") ;
	do {  
		len = getline(&line,&n,efiles) ;
		readEFile(line) ;
/*		printf("%d %d _ %s",len, n,line) ; */
		return ;
	} while (len > 0) ;
	
}
int main(int ac , char **av)
{
	getEList("/data/bergth/4D/testing/sk1/1998/apr/01");
}
