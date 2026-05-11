# Print shell input lines to standard error
set -v
#
# Compile Fortran codes using Intel ifx compiler
ifx -O3 -autodouble -ocnv_nez cnv_nez.for
ifx -O3 -autodouble -ogr_nez gr_nez.for
ifx -O3 -oprepare prepare.f90
ifx -O3 -oresort resort.f90
#
# Prepare AXITRA calculations
./prepare
rm -fr dat
mkdir dat
#
# Compile Fortran codes using older Intel ifort compiler 
#ifort -O3 -autodouble -ocnv_nez cnv_nez.for
#ifort -O3 -autodouble -ogr_nez gr_nez.for
#ifort -O3 -oprepare prepare.f90
#ifort -O3 -oresort resort.f90
