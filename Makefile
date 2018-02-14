
CFLAGS=-g

top : editeve
	editeve  -l -21 -L -20.0

median : median.o
	cc -g  -o median median.o -lm

vf : velfit
	velfit
F= velfit.c ray.o golubc.o 
velfit : $F 
	cc -o velfit -DTEST  $F -lm
	
azimuth.o distance.o editeve.o inverf.o locate.o locate_T.o median.o minvel.o phases.o ray.o rayplot.o readeve.o reltest.o slope.o stations.o travelt.o velfit.o : ray.h

rayplott : rayplot
	rayplot
rayplot : rayplot.o ray.c
	cc -g -DPLOTRAY -o rayplot ray.c rayplot.o -lproj -lm

tt : travelt 
	./travelt -d 3 -x 15

tp : phases 
	./phases -r

tr : 
	cc -g -DTEST readdata.c
	a.out

P =  ray.o readdata.o distance.o stations.o
phases : $P phases.c
	cc -DTEST phases.c $P -o phases -lproj -lm

tl : locate
	locate

R =  ray.o readdata.o distance.o stations.o phases.o reltest.o

trt : reltest
	reltest -r ../geysir/reloc2
reltest : $R
	cc -g -o reltest $R -lproj -lm

L = $P phases.o golubc.o azimuth.o proj.o velfit.o
locate : $L locate.c
	cc -g -o locate  -DTEST $L locate.c -lproj -lm

M =  minvel.o locate.o $L
minvel : $M 
	cc  -g $M -o minvel -lproj -lm

T = travelt.o ray.o

E = readeve.o ray.o  stations.o
readeve : $E
	cc -g $E -lm -o readeve

V = editeve.o locate.o shuffle.o ray.o distance.o stations.o velfit.o golubc.o inverf.o 
$E $V : ray.h Makefile
editeve : $V
	cc -g $V -lm -o editeve -lproj

travelt : $T
	cc ${CFLAGS}  -o travelt $T -lm

t : ray
	./ray

ray : ray.c  ray.h
	cc -g -DDEBUG -o ray ray.c -lm

t22 : travelt
	travelt -n 0 -l 9 -d 7.3 -x 22.1 > jt22  2>&1

t14 : travelt
	travelt -n 0 -l 9  -x 14.6 -d 5.6 > jtest  2>&1

I= rayplot travelt

install : $I
	cp $I ../../bin

clean : 
	rm *.o $I


