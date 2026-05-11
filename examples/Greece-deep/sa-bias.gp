# Spectral acceleration (SA) modeling bias plotted 
# for individual components as a function of a period
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Horizontal layout of plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
set term postscript color enhanced
set output 'sa-bias.ps'
set logscale x
set yrange [-3:3]
set xrange [0.1:10]
set ytics -3,1,3
set grid
unset key
#set ylabel 'ln(SA_{SYN}/SA_{OBS})'
#set tmargin 0.5
#set bmargin 0.5
#set lmargin 5.5
#set rmargin 2.5
NR=18 

#graphLabel="at 1,1 SA BIAS"
#set multiplot layout 1,3 title 'SA BIAS' scale 1,0.4 
set multiplot layout 1,3 scale 1,0.4 
#set multiplot layout 1,3 title 'SA BIAS' 
set arrow 1 from 0.1,0 to 10,0 front nohead lc -1

set style line 1 lt 1 lc rgb "gray60" lw 1 pt 1
set style line 2 lt 1 lc rgb "red"    lw 5 pt 1
set style line 3 lt 2 lc rgb "red"    lw 3 pt 1

# N component
# ~~~~~~~~~~~
set title 'N' offset -8,-2.3
set size 0.35, 0.35
set xlabel 'Period (s)'
set ylabel 'ln(SA_{SYN}/SA_{OBS})'
p for [i=2:NR+1] 'rnkomp.dat' u (1/$1):(column(i)) w l ls 1,\
  '' u (1/$1):(column(NR+2)) w l ls 2,\
  '' u (1/$1):(column(NR+2)-column(NR+3)) w l ls 3 dt 2,\
  '' u (1/$1):(column(NR+2)+column(NR+3)) w l ls 3 dt 2

# E component
# ~~~~~~~~~~~
set title 'E' offset -8,-2.3
set size 0.35, 0.35
set xlabel 'Period (s)'
set ylabel ' { }'
set label 'SA BIAS' at 0.5, 4.5
#unset ylabel
p for [i=2:NR+1] 'rekomp.dat' u (1/$1):(column(i)) w l ls 1,\
  '' u (1/$1):(column(NR+2)) w l ls 2,\
  '' u (1/$1):(column(NR+2)-column(NR+3)) w l ls 3 dt 2,\
  '' u (1/$1):(column(NR+2)+column(NR+3)) w l ls 3 dt 2

# Z component
# ~~~~~~~~~~~
set title 'Z' offset -8,-2.3
set size 0.35, 0.35
set xlabel 'Period (s)'
set ylabel ' { }'
#unset ylabel
unset label
p for [i=2:NR+1] 'rzkomp.dat' u (1/$1):(column(i)) w l ls 1,\
  '' u (1/$1):(column(NR+2)) w l ls 2,\
  '' u (1/$1):(column(NR+2)-column(NR+3)) w l ls 3 dt 2,\
  '' u (1/$1):(column(NR+2)+column(NR+3)) w l ls 3 dt 2
unset multiplot

system "convert -density 300 -rotate 90 -trim sa-bias.ps sa-bias.png"

