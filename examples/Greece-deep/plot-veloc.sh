#!/bin/bash
SCEN=1
Tshift=0.
NP=16384
NPp=4096
#NPp=100
staname=`awk "{if(NR==$1) print \\$4}" stations.dat`
echo $staname
gnuplot << END
set term postscript enhanced color solid 12
set output "plotseis/veloc-$staname.ps"
set multiplot
set size .5,.3
STA1=$1+1
STA2=$1+1
Tshift=$Tshift

set format y "%G"
set ylabel 'Velocity (m/s)'
set xrange [30:120]
#set yrange [-1.e-5:1.e-5]
set origin 0,.7
#unset xtics
set label 'N' offset 23,-2.5 at screen 0.20,1.
plot 'observed/vseisn.dat' u 1:(column(STA1)) every ::(($SCEN-1)*$NP)::($SCEN*$NP-1) notitle w l lt -1 lw .5 lc -1,\
     'vseisn.dat' u (column(1)+Tshift):(column(STA2)) every ::(($SCEN-1)*$NP)::($SCEN*$NP-1) notitle w l lt 1 lc 7
set origin 0,.4
set label 'E' offset 23,-15 at screen 0.20,1.
plot 'observed/vseise.dat' u 1:(column(STA1)) every ::(($SCEN-1)*$NP)::($SCEN*$NP-1) notitle w l lt -1 lw .5 lc -1,\
     'vseise.dat' u (column(1)+Tshift):(column(STA2)) every ::(($SCEN-1)*$NP)::($SCEN*$NP-1) notitle w l lt 1 lc 7
set origin 0,.1
set xlabel 'Time (s)'
set label 'Z' offset 23,-27.5 at screen 0.20,1.
plot 'observed/vseisz.dat' u 1:(column(STA1)) every ::(($SCEN-1)*$NP)::($SCEN*$NP-1) notitle w l lt -1 lw .5 lc -1,\
     'vseisz.dat' u (column(1)+Tshift):(column(STA2)) every ::(($SCEN-1)*$NP)::($SCEN*$NP-1) notitle w l lt 1 lc 7
unset xlabel

set format y "10^{%T}"
set ylabel 'Fourier spectrum (m)'
set xrange [0.01:15]
set yrange [1.e-6:]
set logscale xy
set origin .5,.7
set label 'N' offset 73,-2.5 at screen 0.20,1.
plot 'observed/rvsseisn.dat' u 1:STA1 every ::(($SCEN-1)*$NPp)::($SCEN*$NPp-1) notitle w l lt -1 lw .5 lc -1,\
     'rvsseisn.dat' u 1:STA2 every ::(($SCEN-1)*$NPp)::($SCEN*$NPp-1) notitle w l lt 1 lc 7
set origin .5,.4
set label 'E' offset 73,-15 at screen 0.20,1.
plot 'observed/rvsseise.dat' u 1:STA1 every ::(($SCEN-1)*$NPp)::($SCEN*$NPp-1) notitle w l lt -1 lw .5 lc -1,\
     'rvsseise.dat' u 1:STA2 every ::(($SCEN-1)*$NPp)::($SCEN*$NPp-1) notitle w l lt 1 lc 7
set origin .5,.1
set label 'Z' offset 73,-27.5 at screen 0.20,1.
set xlabel 'Frequency (Hz)'
plot 'observed/rvsseisz.dat' u 1:STA1 every ::(($SCEN-1)*$NPp)::($SCEN*$NPp-1) notitle w l lt -1 lw .5 lc -1,\
     'rvsseisz.dat' u 1:STA2 every ::(($SCEN-1)*$NPp)::($SCEN*$NPp-1) notitle w l lt 1 lc 7
unset xlabel

unset multiplot
END

convert -density 600 -rotate 90 -trim "plotseis/veloc-$staname.ps" "plotseis/veloc-$staname.png"
