#
set -v
# Compile Fortran codes using Intel ifx compiler
ifx -qopenmp -oKK KingKong2.f90 four1.for splie2.for splin2.for spline.for splint.for rlft3.f90 KKgen.f90 JBdistance.f90 filters.for
#ifx -qopenmp -autodouble -oanalyze analyze.f90
ifx -CA -qopenmp -autodouble -oanalyze analyze.f90
ifx -osa-bias-prep sa-bias-prep.f90 
#
# Compile Fortran codes using older Intel ifort compiler 
#ifort -qopenmp -oKK KingKong2.f90 four1.for splie2.for splin2.for spline.for splint.for rlft3.f90 KKgen.f90 JBdistance.f90 filters.for
#ifort -qopenmp -autodouble -oanalyze analyze.f90
#ifort -CA -qopenmp -autodouble -oanalyze analyze.f90
