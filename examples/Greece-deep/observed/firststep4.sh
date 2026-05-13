set -v
#
cp -pv ../../../src/dw/input.dat .
#
# Copy stations.JBdist.dat to current directory
cp -pv ../../../examples/Greece-deep/stations.JBdist.dat .
#
# Make symbolic link between code analyze
ln -s ../../../src/KK-src/analyze analyze
#
# Compile code processseis
ifx -oprocessseis processseis.f90 filters.for
