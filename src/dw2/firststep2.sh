# Print shell input lines to standard error
set -v
#
# Compile Fortran codes using Intel ifx compiler
ifx -O -autodouble -ocnv_nez cnv_nez.for
ifx -O -autodouble -ogr_nez gr_nez.for
ifx -O -oprepare2 prepare2.f90 KKgen.f90 rlft3.f90
ifx -O -oresort2 resort2.f90

# Prepare AXITRA calculations
./prepare2
rm -fr dat
mkdir dat

#ifort -O -autodouble -ocnv_nez cnv_nez.for
#ifort -O -autodouble -ogr_nez gr_nez.for
#ifort -O -oprepare2 prepare2.f90 KKgen.f90 rlft3.f90
#ifort -O -oresort2 resort2.f90
