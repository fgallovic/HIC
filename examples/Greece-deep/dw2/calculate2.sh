# Print shell input lines to standard error
set -v
#
# Number of parallel calculations (generally number of cores)
processors=12
#
# Calculation of Green function
cat sources2.dat | xargs -P$processors -n7 ./gr_nez
#
# Convolution of Green function with moment tensor
cat sources2.dat | xargs -P$processors -n7 ./cnv_nez
#
# Green’s functions are resorted into the NEZsor2.dat
./resort2
#
echo "calculate2.sh Done!"
