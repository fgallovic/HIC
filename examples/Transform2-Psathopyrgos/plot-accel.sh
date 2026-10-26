#!/bin/bash
Tshift=0.
staname=`awk "{if(NR==$1) print \\$4}" stations.dat`
echo $staname
gnuplot << END
set term postscript enhanced color solid 12
set output "plotseis/accel-$staname.ps"
set multiplot
set size .5,.3
STA1=$1+1
STA2=$1+1
Tshift=$Tshift

set format y "%G"
set ylabel 'Acceleration (m/s^2)'
set xrange [0:60]
#set yrange [-1.e-4:]
set origin 0,.7
#unset xtics
set label 'N' offset 23,-2.5 at screen 0.20,1.
plot 'raseisn.dat' u (column(1)+Tshift):(column(STA2)) notitle w l lt 1 lc 7
set origin 0,.4
set label 'E' offset 23,-15 at screen 0.20,1.
plot 'raseise.dat' u (column(1)+Tshift):(column(STA2)) notitle w l lt 1 lc 7
set origin 0,.1
set xlabel 'Time (s)'
set label 'Z' offset 23,-27.5 at screen 0.20,1.
plot 'raseisz.dat' u (column(1)+Tshift):(column(STA2)) notitle w l lt 1 lc 7
unset xlabel

set ylabel 'Fourier spectrum (m/s^2)'
set xrange [0.01:15]
set yrange [1.e-6:]
set logscale xy
set origin .5,.7
set label 'N' offset 73,-2.5 at screen 0.20,1.
plot 'rasseisn.dat' u 1:STA2 notitle w l lt 1 lc 7
set origin .5,.4
set label 'E' offset 73,-15 at screen 0.20,1.
plot 'rasseise.dat' u 1:STA2 notitle w l lt 1 lc 7
set origin .5,.1
set label 'Z' offset 73,-27.5 at screen 0.20,1.
set xlabel 'Frequency (Hz)'
plot 'rasseisz.dat' u 1:STA2 notitle w l lt 1 lc 7
unset xlabel

unset multiplot
END

convert -density 600 -rotate 90 -trim "plotseis/accel-$staname.ps" "plotseis/accel-$staname.png"
