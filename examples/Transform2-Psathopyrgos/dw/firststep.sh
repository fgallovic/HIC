set -v
#
# Make symbolic links between executable files
ln -s ../../../src/prepare prepare
ln -s ../../../src/gr_nez gr_nez
ln -s ../../../src/cnv_nez cnv_nez
ln -s ../../../src/resort resort
#
# Prepare AXITRA calculations
./prepare
rm -fr dat
mkdir dat
