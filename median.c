#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ray.h"
char *selfdoc=" median file1 file2 ..... filen\n\
\n\
\n\
\n\
\n\
\n\
\n" ;

#define NCOL 5
double data[10000] ;
char line[1000] ;
int nData,nLines,nFiles, nLevel ;

int compareDouble( const void *p1, const void *p2 )
{
        double d1,d2 ;
        d1 = *( double *) p1 ;
        d2 = *( double *) p2 ;
        if( d2 < d1 ) return 1 ;
        if( d2 > d1 ) return -1 ;
        return 0 ;
}
void printData() 
{
	int n,j,i ;
	n = nLevel*NCOL*nFiles ;
	for( i = 0 ; i < n ; i++) {
		j = i % NCOL ;
		if( j == 0 ) printf("\n") ;
		printf("%9.3f",data[i]) ;
	}
	printf("\n") ;
}
void doMedian()
{
	double work[1000],median ;
	int stride, offset ;
	int iLevel, iFile, iCol  ;
	stride = NCOL*nLevel ;
	for( iLevel = 0 ; iLevel < nLevel ; iLevel++) {
		printf("RV" ) ;
		for( iCol = 0 ; iCol < NCOL ; iCol++ ) {
			offset = iLevel*NCOL + iCol ;
			for(iFile = 0 ; iFile < nFiles ; iFile++)
				work[iFile] = data[offset + iFile*stride] ;
			qsort(work,nFiles,sizeof(double),compareDouble) ;
			median = 0.5*(work[nFiles/2] + work[(nFiles+1)/2]) ;
			printf("%9.3f",median) ;
		}
		printf("\n") ;
	}
}
int main(int ac, char **av)
{
	int i,iLine,j ;
	char *fn ;
	FILE *file ;
	char *res ;
	double *dp ;
	nFiles = ac - 1 ;
	if( nFiles < 2 ) {
		printf(selfdoc) ;
		exit(0) ;
	}	
	dp = data ;
	for( i = 1 ; i < ac ; i++) {
		iLine = 0 ;
		fn = av[i] ;
		file = fopen(fn,"r") ;
		if( NULL == file) {
			printf("Error opening %s\n",fn) ;
			exit(-1) ;
		}
		while (res = fgets(line,1000,file)) 
		    if( (*res == 'R') && ( res[1] == 'V')) {
			for( j = 0 ; j < NCOL ; j++ ) *dp++ = atof(line+4+9*j) ;
			iLine++ ;
#ifdef DEBUG
			printf("%d %s %s",iLine,fn,line ) ;
#endif
		}
	}
	nData = dp-data ;
	nLines = nData/NCOL ;
	nLevel = nLines/nFiles ;
#ifdef DEBUG
	printf("%d %d\n",nData,nLevel*nFiles*NCOL ) ;
	printData() ;
#endif
	doMedian() ;
}
