# PGD, PGV and PGA from horizontal components as a function of azimuth and distance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set term postscript color enhanced 10
set output 'pg-d-v-a.ps'

set grid
unset key
set logscale y
set style circle radius 1

set multiplot layout 3,2

# Peak ground displacements (PGD)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set ylabel 'PGD (m)' offset 1.5,0
set xrange [-180:180]
p 'rpgd.dat' u ($10*180/pi):8 w p pt 7 lc 7 title 'Synt'

unset ylabel
unset xlabel
set xrange [1:100]
set logscale x
p 'rpgd.dat' u 11:8 w p pt 7 lc 7 title 'Synt'
unset logscale x

# Peak ground velocities (PGV)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set ylabel 'PGV (m/s)' offset 1.5,0
set xrange [-180:180]
p 'rpgv.dat' u ($10*180/pi):8 w p pt 7 lc 7 title 'Synt'

set xrange [1:100]
set logscale x
unset ylabel
unset xlabel
p 'rpgv.dat' u 11:8 w p pt 7 lc 7 title 'Synt',\
  'PGV-Bea-21.txt' u 1:($2) w l lc 1 lw 3 title 'Bea2021',\
  'PGV-Bea-21.txt' u 1:($3) w l lc 1 lw 3 dt 2 notitle,\
  'PGV-Bea-21.txt' u 1:($4) w l lc 1 lw 3 dt 2 notitle
unset logscale x

# Peak ground accelerations (PGA)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set ylabel 'PGA (m/s^2)'
set xlabel 'Azimuth (deg)'
set xrange [-180:180]
p 'rpga.dat' u ($10*180/pi):8 w p pt 7 lc 7 title 'Synt'

unset ylabel
set key b l inside
set xlabel 'J-B distance (km)'
set xrange [1:100]
set logscale x
p 'rpga.dat' u 11:8 w p pt 7 lc 7 title 'Synt',\
  'PGA-Bea-21.txt' u 1:($2) w l lc 1 lw 3 title 'Bea2021',\
  'PGA-Bea-21.txt' u 1:($3) w l lc 1 lw 3 dt 2 notitle,\
  'PGA-Bea-21.txt' u 1:($4) w l lc 1 lw 3 dt 2 notitle
unset multiplot

system "convert -density 300 -rotate 90 -trim pg-d-v-a.ps pg-d-v-a.png"
