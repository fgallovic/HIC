# PGD, PGV and PGA from horizontal components as a function of azimuth and distance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set term postscript color enhanced 10
set output 'pg-d-v-a.ps'

#set size 1.,0.33
set grid
unset key
#set arrow 1 from -155, 0.01 to -155,1. nohead 
#set arrow 2 from 60, 0.01 to 60,1. nohead 
set logscale y
set yrange [0.0001:]
#set rmargin 2.
#set lmargin 3.
#set tmargin 1.
#set bmargin 1.
set style circle radius 1
depth=88.6

set multiplot layout 3,2

# Peak ground displacements (PGD)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set ylabel 'PGD (m)' offset 1.5,0
set xrange [-180:180]
#set arrow 1 from -155, 0.01 to -155,10. nohead 
#set arrow 2 from 60, 0.01 to 60,10. nohead 
p 'rpgd.dat' u ($10*180/pi):8 w p pt 7 lc 7 title 'Synt',\
  'observed/rpgd.dat' u ($10*180/pi):8 w p pt 7 ps 0.75 lc -1 title 'Obs'
#unset arrow 1
#unset arrow 2

#set ylabel 'PGD (m)'
unset ylabel
unset xlabel
#set xlabel 'Distance (km)'
set xrange [80:400]
set xtics ("100" 100,"200" 200,"300" 300) 
set logscale x
p 'rpgd.dat' u 11:8 w p pt 7 lc 7 title 'Synt',\
  'observed/rpgd.dat' u 11:8 w p pt 7 ps 0.75 lc -1 title 'Obs'
unset logscale x

# Peak ground velocities (PGV)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set ylabel 'PGV (m/s)' offset 1.5,0
#set xlabel 'Azimuth (deg)'
set xrange [-180:180]
#set arrow 1 from -155, 0.01 to -155,10. nohead 
#set arrow 2 from 60, 0.01 to 60,10. nohead 
#set title 'PGV vs Angle'
p 'rpgv.dat' u ($10*180/pi):8 w p pt 7 lc 7 title 'Synt',\
  'observed/rpgv.dat' u ($10*180/pi):8 w p pt 7 ps 0.75 lc -1 title 'Obs'
#unset arrow 1
#unset arrow 2

set xrange [80:400]
set xtics ("100" 100,"200" 200,"300" 300) 
set logscale x
unset ylabel
unset xlabel
p 'rpgv.dat' u  (sqrt($2**2+$3**2+depth**2)):($8*100) w p pt 7 lc 7 title 'Synt',\
  'observed/rpgv.dat' u  (sqrt($2**2+$3**2+depth**2)):($8*100) w p pt 7 ps 0.75 lc -1 title 'Obs',\
  'pgv-back-arc-avg.out' u 1:($2) w l lc rgb "dark-spring-green" lw 3 title 'backarc GMM',\
  'pgv-back-arc-avg.out' u 1:($3) w l lc rgb "dark-spring-green" lw 3 dt 2 notitle,\
  'pgv-back-arc-avg.out' u 1:($4) w l lc rgb "dark-spring-green" lw 3 dt 2 notitle,\
  'pgv-along-arc-avg.out' u 1:($2) w l lc rgb "#D55E00" lw 3 title 'forearc GMM',\
  'pgv-along-arc-avg.out' u 1:($3) w l lc rgb "#D55E00" lw 3 dt 2 notitle,\
  'pgv-along-arc-avg.out' u 1:($4) w l lc rgb "#D55E00" lw 3 dt 2 notitle
unset logscale x

# Peak ground accelerations (PGA)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set ylabel 'PGA (m/s^2)'
set xlabel 'Azimuth (deg)'
set xrange [-180:180]
set yrange [0.001:]
#set arrow 1 from -155, 0.01 to -155,100. nohead 
#set arrow 2 from 60, 0.01 to 60,100. nohead 
p 'rpga.dat' u ($10*180/pi):8 w p pt 7 lc 7 title 'Synt',\
  'observed/rpga.dat' u ($10*180/pi):8 w p pt 7 ps 0.75 lc -1 title 'Obs'
#unset arrow 1
#unset arrow 2

unset ylabel
set key b l inside
#set ylabel 'PGA (m/s^2)'
set xlabel 'Distance (km)'
set xrange [80:400]
set xtics ("100" 100,"200" 200,"300" 300) 
set logscale x
#set offset 10,10,10,10
p 'rpga.dat' u  (sqrt($2**2+$3**2+depth**2)):($8*100) w p pt 7 lc 7 title 'Synt',\
  'observed/rpga.dat' u  (sqrt($2**2+$3**2+depth**2)):($8*100) w p pt 7 ps 0.75 lc -1 title 'Obs',\
  'pga-back-arc-avg.out' u 1:($2) w l lc rgb "dark-spring-green" lw 3 title 'backarc GMM',\
  'pga-back-arc-avg.out' u 1:($3) w l lc rgb "dark-spring-green" lw 3 dt 2 notitle,\
  'pga-back-arc-avg.out' u 1:($4) w l lc rgb "dark-spring-green" lw 3 dt 2 notitle,\
  'pga-along-arc-avg.out' u 1:($2) w l lc rgb "#D55E00" lw 3 title 'forearc GMM',\
  'pga-along-arc-avg.out' u 1:($3) w l lc rgb "#D55E00" lw 3 dt 2 notitle,\
  'pga-along-arc-avg.out' u 1:($4) w l lc rgb "#D55E00" lw 3 dt 2 notitle
unset multiplot

system "convert -density 300 -rotate 90 -trim pg-d-v-a.ps pg-d-v-a.png"
