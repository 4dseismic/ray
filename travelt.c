/* 
	Program to print out traveltime data for 1d linear gradeint layered models.
	Einar Kjartransson, 2016 
*/
char *selfDoc = "\
travelt options\n\
Program uses raytracing throuch 1D earth that consists of layers where velocity varies linearly with depth.\n\
	Options:\n\
-f	velfile	   silp.vel	File that contains the velocity function.\n\
-l	logLevel		Determines the amount of information printed.\n\
-d	sourceDepth	4.5	Source depth of earthquake.\n\
-n	nPoint	 	 10	Number of points top plot in each raytracing domain.\n\
-b	bottomDepth	 15	Maximum turning depth for rays.\n\
-v	velReduce	  7	Velocity to use to reduce traveltime. Value of 0 yields full traveltime.\n\
-o	outputFile		Write to outputFile. If not specified name is generated form option valuses.\n\
-D				Raytrace for depth pases.\n\
-r	n dz			Velocity function will be resampled before use, using cubic splines.\n\
-X	xEnd xCount		Use iteration to output values for evently spaces distance values.\n\
";
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define _GNU_SOURCE
#include <fenv.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include "ray.h"

char *velFile = "silp.vel" ;
double sourceDepth = 4.5 ;
double bottomDepth = 15 ;
double velReduce = 7.0 ;
int nPoint = 10 ;
int depthPhaseFlag ;
double xSolve,zSolve ;
char *tableFile ;
double xEnd ;
int xCount ;

int nResample = 0 ;
double dzResample ; 
char outputName[300] ;
VelModel m ;

void doXList() 
{
	int i ;
	double x,t,dxdp,p,dtdx ;
	FILE *fd ;
	fd = fopen("travelt.xlist","w") ;
	for( i = 0 ; i < xCount ; i++) {
		x = (i + 0.5) * xEnd / xCount ;
		t = timeFromDist(&m,x,sourceDepth,&p,&dtdx,&dxdp ) ;
		if( 0.001 >= fabs(dxdp) ) 
		  fprintf(fd,"%10.6f  %10.6f %10.6f %10.6f %10.6f\n",x,t,p,dtdx,dxdp) ;
	}
	fclose(fd) ;
}
void readFromTable()
{
	FILE *table ;
	char line[305] ;
	int i,len ;
	double x,z,t ;
	double tt,dtdx,dxdp,p,dt,sum1,sigma,res,sum2,sigma2,drt ;
	table = fopen(tableFile,"r") ;
	sum1 = 0.0 ;
	sum2 = 0.0 ;
	i = 0 ;
	len = 300 ;
	while(  fgets( line,len,table) ) {

		sscanf(line,"%lf %lf %lf %lf",&x,&z,&t,&res) ;
		tt = timeFromDist(&m,x,z,&p,&dtdx,&dxdp) ;
		drt = tt-t+res ;
		dt = tt-t ;
		printf("%10.4f %10.4f %10.4f %10.4f %8.3f %8.3f %10.5f %s",x,z,t,tt,dt,drt,1.0/dxdp,line) ;
		i++ ;
		sum1 += dt*dt ;
		sum2 += drt*drt ;
	}
	sigma = sqrt(sum1/i ) ;
	sigma2 = sqrt(sum2/i ) ;
	printf("n=%d sigma=%10.4f sigma2=%8.3f\n",i,sigma,sigma2) ;
}
void doIt() 
{
	double pMax, pBottom, p,pp , dp, x,t, reduce ;
	double dtdx,dpdx,p1,x1,dxdp ;
	int il,i ;
	Mode mmode ;
	FILE *ofd ;
	ofd = fopen(outputName,"w") ;
	if( velReduce > 0.0 ) reduce = 1.0/velReduce ; else reduce = 0.0 ; 
	initVelModel(200, &m  ) ;
	readVelModel(velFile, &m) ;
	if( nResample > 1 ) { m = resampleVelModel(&m,dzResample,nResample); }  
	if( shLogLevel >= 3 ) printVelModel( &m) ;
	pMax = 1.0/velZ(sourceDepth,&m,&il ) ;
	mmode = RayDown ;
	if( depthPhaseFlag ) mmode = DepthPhase ;
	if( sourceDepth <= 0.0 ) mmode = Surface ; 
	if ( mmode < Surface)  for( i = 0 ; i < nPoint ; i++) {
		p1 = (i + 0.8) *pMax / nPoint ;
		x1 = traceUD( RayUp, p1, sourceDepth, &m, &t) ;
		p = (i + 0.81) *pMax / nPoint ;
		x = traceUD( RayUp, p, sourceDepth, &m, &t) ;
		dpdx = ( p1-p) / (x1-x) ;
		fprintf(ofd,"%10.5f %10.5f %6d %10.6f %10.6f\n",x,t-x*reduce,i,p,dpdx) ;
	}
	pBottom = 1.0/velZ(bottomDepth,&m,&il) ;
	dp = pMax - pBottom ;
	for( i = 0 ; i < nPoint ; i++) {
/*		p = pMax - (i+0.1)*dp/nPoint ;  */
		p1 = 1.0/velZ( sourceDepth + (i+0.21)*(bottomDepth-sourceDepth)/nPoint,&m,&il);   
		x1 = traceUD( mmode , p1, sourceDepth, &m, &t) ; 
		p = 1.0/velZ( sourceDepth + (i+0.2)*(bottomDepth-sourceDepth)/nPoint,&m,&il);   
		x = traceUD( mmode , p, sourceDepth, &m, &t) ; 
		dpdx = ( p1-p) / (x1-x) ;
		fprintf(ofd,"%10.5f %10.5f %6d %10.6f %10.6f %10.4f\n",x,t-x*reduce,i,p,dpdx,zBottom) ;
	} 
	if( 0.0 != xSolve) t = timeFromDist(&m,xSolve,sourceDepth,&pp,&dtdx,&dxdp) ;
	if( tableFile) readFromTable() ;
	fclose(ofd) ;
}
void makeOutputName( int i , char *a ) 
{
	char  *b ;
	b = outputName ; 
	if(*a == '.') a++ ;
	if(*a == '/') a++ ;
	b=strcpy(b,"P") ; b += 1 ;
	while ( i > 0 ) {
		*b = *a++ ;
		if( *b++ == 0 ) { b-- ; i-- ; }
	} 
}
int main(int ac , char **av ) 
{
	extern char *optarg ;
	extern int optind ;
	int cc ;
	feenableexcept(FE_INVALID) ;
	while( EOF != (cc = getopt(ac,av,"f:l:d:n:b:v:o:sr:x:t:X:DhH?"))) {
		switch(cc) {
		case 'f':	velFile=optarg ; break ;
		case 'l' :	shLogLevel = atoi(optarg) ; break ;
		case 'd' :      sourceDepth = atof(optarg) ; break ;
		case 'n' : 	nPoint = atoi(optarg) ; break ;
		case 'b' : 	bottomDepth = atof(optarg) ; break ;
		case 'v' :	velReduce = atof(optarg) ; break ;
		case 'o' :	strcpy( outputName,optarg) ; break ;
		case 's' : 	splineTest() ; break ;
		case 'r' :	nResample = atoi(optarg) ; dzResample = atof(av[optind++]) ; break ;
		case 'x' :      xSolve = atof(optarg) ; break ;
		case 't' :	tableFile = optarg ; break ;
		case 'D' :	depthPhaseFlag = 1 ; break ;
		case 'X' : 	xEnd = atof(optarg) ; xCount =  atof(av[optind++]) ; break ;
		default : printf(selfDoc) ; break ;
	}}
	if( sourceDepth <= 0.0 ) depthPhaseFlag = 0 ;
	if( 0 == *outputName )  makeOutputName(ac,*av) ; 
	fprintf(stderr,"_%s_\n",outputName) ;
	doIt() ;
	if( xCount ) doXList() ;
	return 0 ; 
}
