set -v
#
# Copy files to current directory
cp -pv dw2/input.dat .
cp -pv dw2/stations.dat .
cp -pv dw2/fault.dat .
#
# Make symbolic links between files
ln -s ../../src/KK KK
ln -s ../../src/analyze analyze
ln -s ../../src/sa-bias-prep sa-bias-prep
ln -s dw/XYGreen.dat XYGreen.dat
ln -s dw/NEZsor.dat NEZsor.dat
ln -s dw2/NEZsor2.dat NEZsor2.dat

