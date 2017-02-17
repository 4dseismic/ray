#!/bin/bash

d=6.4
n=`basename $0 .sh`
o=$n.svg
t=`echo $n | tr '_' ' ' `
v=silp.vel

cat <<EOF > gplot.input
set terminal svg size 1400,600 font ",24"
set output "$o"
set grid
set yrange[0:3]
set xrange [0:115]
#set size ratio -1
set title "P, S and P depth phase. Source at $d km"
set xlabel "distance (km)"
set ylabel "reduced travel time (s)"
#set y2tics 0 0.25
set xtics 5
#set y2label "reduced travel time (s)"
#set link y2 via y/8.0 inverse y*8.0
set font ",20"
unset key
plot \\
EOF
l="with lines"

awk ' {  print $2,$1 } ' dump.vel > velocity
#cp dump.vel velocity

./travelt -o jj -f $v -r 370 0.05 -d $d -n 1300 -v 6.5 -b $i 
echo >> jj
./travelt -o jjj -f sils.vel -r 370 0.05 -d $d -n 1300 -v 6.5 -b $i 
cat jjj >> jj
echo >> jj
./travelt -o jjj -f $v -r 370 0.05 -D -d $d -n 1300 -v 6.5 -b $i 
cat jjj >> jj

echo \"jj\"  $l , \\ >> gplot.input
gnuplot < gplot.input

display $o
