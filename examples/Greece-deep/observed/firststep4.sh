set -v
#
cp -pv ../synthetic/dw/input.dat .
#
# Copy stations.JBdist.dat to current directory
cp -pv ../synthetic/KK/stations.JBdist.dat .
#
# Make symbolic link between code analyze
ln -s ../synthetic/KK-src/analyze analyze
#
# Compile code processseis
ifx -oprocessseis processseis.f90 filters.for
