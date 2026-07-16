# Print shell input lines to standard error
set -v
#
# Compile Fortran codes using Intel ifx compiler
#
# dw + dw2
# ^^^^^^^^
ifx -O3 -autodouble -o cnv_nez cnv_nez.for
ifx -O3 -autodouble -o gr_nez gr_nez.for
ifx -O3 -o prepare prepare.f90
ifx -O3 -o resort resort.f90
# dw2
# ^^^
ifx -O -o prepare2 prepare2.f90 KKgen.f90 rlft3.f90
ifx -O -o resort2 resort2.f90
# KK 
# ^^
ifx -qopenmp -o KK KingKong2.f90 four1.for splie2.for splin2.for spline.for splint.for rlft3.f90 KKgen.f90 JBdistance.f90 filters.for
#ifx -qopenmp -autodouble -o analyze analyze.f90
ifx -CA -qopenmp -autodouble -o analyze analyze.f90
ifx -o sa-bias-prep sa-bias-prep.f90 
# observed
# ^^^^^^^^
ifx -o processseis processseis.f90 filters.for
