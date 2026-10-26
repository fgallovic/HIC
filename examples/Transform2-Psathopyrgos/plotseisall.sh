#!/bin/bash
rm -rf plotseis
mkdir plotseis
processors=12
stations=42
seq $stations | xargs -P$processors -n1 ./plot-displ.sh
seq $stations | xargs -P$processors -n1 ./plot-veloc.sh
seq $stations | xargs -P$processors -n1 ./plot-accel.sh
#
#seq $stations | xargs -P$processors -n1 ./plot-displ-synt.sh
#seq $stations | xargs -P$processors -n1 ./plot-veloc-synt.sh
#seq $stations | xargs -P$processors -n1 ./plot-accel-synt.sh
echo "Done"
