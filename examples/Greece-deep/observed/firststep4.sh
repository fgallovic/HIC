set -v
#
# Copy files to current directory
cp -pv ../dw/input.dat .
cp -pv ../analyze.in .
cp -pv ../frequencies.txt .
cp -pv ../stations.JBdist.dat .
#
# Make symbolic links between codes
ln -s ../../../src/analyze analyze
ln -s ../../../src/processseis processseis
