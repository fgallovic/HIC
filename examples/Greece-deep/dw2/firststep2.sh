set -v
#
# Make symbolic links between executable files
ln -s ../../../src/prepare2 prepare2
ln -s ../../../src/gr_nez gr_nez
ln -s ../../../src/cnv_nez cnv_nez
ln -s ../../../src/resort2 resort2
#
# Prepare AXITRA calculations
./prepare2
rm -fr dat
mkdir dat
