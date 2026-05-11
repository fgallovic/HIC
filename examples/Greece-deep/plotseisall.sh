#!/bin/bash
rm -rf plotseis
mkdir plotseis
processors=12
stations=18
seq $stations | xargs -P$processors -n1 ./plot-displ.sh
seq $stations | xargs -P$processors -n1 ./plot-veloc.sh
seq $stations | xargs -P$processors -n1 ./plot-accel.sh
echo "Done"
