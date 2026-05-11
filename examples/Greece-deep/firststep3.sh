set -v
#
cp -pv ../../src/dw2/input.dat .
cp -pv ../../src/dw2/stations.dat .
cp -pv ../../src/dw2/fault.dat .
#cp -pv ../../src/dw2/GreeceDeep-slip2D.dat .

cd ../../src/KK-src
./compile.sh
cd -
#
# Make symbolic links between files
ln -s ../../src/KK-src/KK KK
ln -s ../../src/KK-src/analyze analyze
ln -s ../../src/KK-src/sa-bias-prep sa-bias-prep
ln -s ../../src/dw/XYGreen.dat XYGreen.dat
ln -s ../../src/dw/NEZsor.dat NEZsor.dat
ln -s ../../src/dw2/NEZsor2.dat NEZsor2.dat

