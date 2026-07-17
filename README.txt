Hybrid integral-composite (HIC) broadband kinematic source modeling
===================================================================

Hybrid Integral-Composite (HIC) technique is used for the broadband 
ground motion simulations (Gallovic and Brokesova, 2007). 
The method combines low-frequency deterministic and high-frequency 
stochastic modeling up to frequencies of engineering interest.

Overview
********
This manual describes how to use open-source software package HIC.
The present version of HIC is suitable for UNIX systems, and we test
and use the software package under Linux. For Windows users, modification 
of Linux shell scripts needs to be performed. Fortran and Gnuplot  
commands should be operating system independent. We tested all 
calculations on Ubuntu 22.04 using Intel Fortran compiler ifx 
(replacing the older ifort).
   Broadband ground motion simulations performed by HIC modelling codes
are demonstrated on an example of an Mw5.8, 29.8. 2014, intermediate-depth 
earthquake in the Hellenic slab. 
   At first we calculate synthetic data for 18 selected seismic stations
and then we synthesize results (directory 'examples/Greece-deep' and
subdirectories 'dw' and 'dw2'). Afterwards we prepare observed velocity 
seismograms from 18 selected real seismic stations for comparison with 
synthetic data (directory 'examples/Greece-deep/observed'). Then we use 
several scripts to plot the results (directory 'examples/Greece-deep').

We shall describe directory structures:

Directory 'src'
~~~~~~~~~~~~~~~
Directory contains source codes for HIC calculations and script
for the compilation of source codes.

Directory 'examples/Greece-deep'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Directory where we synthesize results from 'dw' and 'dw2' subdirectories. 
Output files contain displacement, velocity, acceleration waveforms, 
Fourier amplitude spectra, peak ground motion values (displacements, 
velocities, accelerations), Housner spectral intensity.
Directory 'examples/Greece-deep' also contains scripts for plotting the final 
results combined with observed data.

Directory 'examples/Greece-deep' has following subdirectories:

'dw' - (discrete wavenumber) directory where we apply the integral approach 
     at low frequencies. 
     Calculation is based on modified Axitra code based on discrete wavenumber 
     method providing full-wavefield Green's functions in 1D layered media
     (Bouchon, 1981; Cotton and Coutant, 1997).
     The code was originally written by O. Coutant.
     Additional modifications have been made by J. Zahradnik, J. Burjanek 
     and F. Gallovic.

'dw2' - here we use the composite approach at high frequencies. 
      Calculation uses the same Axitra code based on discrete wavenumber 
      method. Some input parameters are different in comparison with
      directory 'dw'.

'observed' - here we prepare observed velocity seismograms from 18 selected 
     real seismic stations for comparison with synthetic data.

How to use the code
*******************

1. Main input data files (common for several directories)
   ~~~~~~~~~~~~~~~~~~~~~
   We shall demonstrate the structure of input files on data prepared
   for event Mw5.8, 29.8. 2014, intermediate-depth earthquake in 
   the Hellenic slab.

   input.dat (in directory 'examples/Greece-deep' and subdirectories 'dw', 'dw2', 'observed') 
   ^^^^^^^^^
   File contains parameters specifying the number of computing frequencies, 
   length of calculated seismograms, number of receivers, and source definition; 
   the latter includes the rectangular source geometry (position, length, width, 
   discretization for Green's functions), as well as moment, strike, dip, rake, 
   rupture velocity. The parameters are read by Fortran executable codes in 
   different directories selectively, i.e., not all parameters are used by
   each executable.

   crustal.dat (in subdirectories 'examples/Greece-deep/dw and dw2') 
   ^^^^^^^^^^^
   Describes a 1D velocity model. We specify the number of horizontal layers, 
   P and S wave velocities, density, and quality factors Qp, Qs, constant in 
   the layers.

   stations.dat (in directory 'examples/Greece-deep' and subdirectories 'dw', 'dw2')
   ^^^^^^^^^^^^
   Contains specification of Cartesian coordinates x, y, z and names of receivers.

   frequencies.txt (directories 'examples/Greece-deep', 'examples/Greece-deep/observed')
   ^^^^^^^^^^^^^^^
   Contains number of frequencies and frequencies of SDF oscillator.
   SDF (Single Degree of Freedom) oscillator represents a simplified model 
   of a structure that vibrates in one dimension. 

   processseis.in (directory 'examples/Greece-deep/observed')
   ^^^^^^^^^^^^^^
   Summarizes various information about observed seismograms: coordinates 
   of seismic stations, origin time, sampling interval, filename, etc.

2. Perform calculations (step by step)
   ~~~~~~~~~~~~~~~~~~~~
   Script 'compile.sh' use commands for Intel Fortran ifx compiler. In the case 
   you need to apply another Fortran compiler modify specific lines. Calculations 
   are set to 12 cores in scripts and calculation times correspond to PC with 
   12 core processor (AMD Ryzen 9 3900X) and SSD NVMe disc. Memory requirements 
   are low in our example.

2.1 Compilation of source codes
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    - directory 'src'
    
    In the directory 'src' run shell script 'compile.sh', that compiles all 
    Fortran codes needed for further calculations.

2.2 Synthetic calculations for real stations
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    - directory 'examples/Greece-deep', subdirectories 'dw', 'dw2'
    
2.2.1 In the directory 'dw' (integral approach at low frequencies) run at first 
      shell script 'firststep.sh', that creates symbolic links to executable files
      compiled in directory 'src' and runs code 'prepare' for preparation of 
      the AXITRA calculations, in particular it prepares list of elementary sources 
      covering the rupture in regular grid. 'firststep.sh' also deletes subdirectory 
      'dat' if exists from previous calculation and makes a new empty subdirectory 
      'dat'.
      
      After successful execution of 'firststep.sh', run shell script 'calculate.sh'
      that starts parallel loop of processes using 'xargs' with option to set 
      the number of processors used for calculation. For each elementary source 
      the codes 'gr_nez' and 'cnv_nez' are run automatically. Intermediate 
      results including Green's functions for the individual elementary sources 
      are stored in the 'dat' directory. Finally, the Green's functions are resorted 
      by 'resort' into the 'NEZsor.dat' file. Note that the order of the Green's 
      functions is such that the outer loop is over stations. Thus if more crustal 
      models are to be considered, the 'NEZsor.dat' files for the individual subsets 
      of stations and then simply appended one after each other.
      Calculation time is 5 minutes.

2.2.2 In the directory 'dw2' (composite approach at high frequencies) run, 
      analogically as in 'dw', at first shell script 'firststep2.sh', that creates 
      symbolic links to executable files in directory 'src' and runs code 
      'prepare2' for preparation of the AXITRA calculations, in particular 
      it prepares list of elementary sources covering the rupture in regular grid.

      After successful execution of 'firststep2.sh', run shell script 'calculate2.sh'
      that starts parallel loop of processes using 'xargs' with option to set 
      the number of processors used for calculation. For each elementary source 
      the codes 'gr_nez' and 'cnv_nez' are run automatically. Intermediate 
      results including Green's functions for the individual elementary sources 
      are stored in the 'dat' directory. Finally, the Green's functions are resorted 
      by 'resort2' into the 'NEZsor.dat' file. Note that the order of the Green's 
      functions is such that the outer loop is over stations. Thus if more crustal 
      models are to be considered, the 'NEZsor.dat' files for the individual subsets 
      of stations and then simply appended one after each other.
      Calculation time is 1 hour and 12 minutes.

2.2.3 Run script 'firststep3.sh' in directory 'examples/Greece-deep'. The script 
      copies necessary input files from directory 'dw2', creates symbolic links 
      to executable codes in directory 'src' and to final output data files 
      in directories 'dw' and 'dw2'. Afterwards, run script 'calculate3.sh'. 
      The script runs code 'KK' to synthesize final results calculated in 'dw' 
      and 'dw2' directories and runs code 'analyze' to produce final set of output 
      files: displacement, velocity, acceleration waveforms, Fourier amplitude 
      spectra, peak ground motion values (displacements, velocities, accelerations), 
      Housner spectral intensity, spectral velocities and accelerations, Arias 
      intensity.  
      Calculation time is few seconds.
      
      'analyze' output files used in plotting: 
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      rdseisn.dat,   rdseise.dat,   rdseisz.dat  ... displacement
      rdsseisn.dat,  rdsseise.dat,  rdsseisz.dat ... displacement spectrum
      vseisn.dat,    vseise.dat,    vseisz.dat   ... velocity
      rvsseisn.dat,  rvsseise.dat,  rvsseisz.dat ... velocity spectrum
      raseisn.dat,   raseise.dat,   raseisz.dat  ... acceleration
      rasseisn.dat,  rasseise.dat,  rasseisz.dat ... acceleration spectrum
      rassseisn.dat, rassseise.dat, rassseisz.dat... smoothed acceleration spectrum
      rpgd.dat, rpgv.dat, rpga.dat ................. peak ground motion values
      rapseisn.dat, rapseise.dat, rapseisz.dat ..... spectral acceleration

2.3 Preparation of observed data for comparison with synthetics
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    - directory 'examples/Greece-deep/observed'

    Go to directory 'examples/Greece-deep/observed'. In our example we use data 
    (velocity seismograms) from 18 stations located in the Hellenic subduction 
    area. At first run script 'firststep4.sh'. The script copies several necessary 
    input files from other directories, one of them is file stations.JBdist.dat 
    (Joyner-Boore distance from a fault) that is prepared during calculations
    in point 2.2.3 in directory 'examples/Greece-deep'. 'firststep4.sh' also
    creates symbolic links to codes 'analyze' and 'processseis' located
    in directory 'src'.

    Script 'calculate4.sh' runs codes 'processseis' and 'analyze'. 
    Code 'processseis' converts measured velocity seismograms to HIC format 
    according to rules specified in data file 'processseis.in'. Execution 
    of code 'analyze' yields similar output files as for synthetic calculations 
    that are used for next processing and plotting.
    Calculation time is few seconds.

2.4 Plotting results
    ~~~~~~~~~~~~~~~~
    - directory 'examples/Greece-deep'

    Go back to the directory 'examples/Greece-deep'. Run script 'plotseisall.sh',
    based on shell and Gnuplot commands, that creates *.ps (PostScript) 
    and *.png (Portable Network Graphics) files with displacement, velocity, 
    acceleration waveforms and appropriate Fourier amplitude spectras for 
    selected 18 seismic stations, both observed and synthetic (in subdirectory 
    'plotseis').
    
    Run code './sa-bias-prep' from command line. The code prepares input data 
    for gnuplot script 'sa-bias.gp' that displays comparison of synthetic 
    and observed data in the form of modeling bias ln(SAsyn/SAobs), where 
    ln is natural logarithm, SAsyn and SAobs are synthetic and observed spectral 
    accelerations, respectively (Pitarka et al., 2021)

    Gnuplot script 'pg-d-v-a.gp' plots peak ground displacements (PGD), 
    velocities (PGV), and accelerations (PGA) from horizontal components 
    as a function of azimuth and distance. We compare PGV and PGA values
    with empirical GMM (Skarlatoudis et al., 2013).

------------------------------------------------------------------------

References

Bouchon, M. (1981). A simple method to calculate Green's functions
for elastic layered media. Bull. Seism. Soc. Am., 71, 959-971.
doi: 10.1785/BSSA0710040959

Cotton, F. & Coutant, O. (1997). Dynamic stress variations due to shear 
faults in a plane-layered medium. Geophys. J. Int., Vol 128, 676-688.
doi:10.1111/j.1365-246X.1997.tb05328.x
 
Gallovic, F. & Brokesova, J. (2007). Hybrid k-squared Source Model 
for Strong Ground Motion Simulations: Introduction. Phys. Earth Planet. 
Interiors, 160, 34-50. doi: 10.1016/j.pepi.2006.09.002.

Pitarka, A. et al. (2021). Refinements to the Graves-Pitarka kinematic 
rupture generator, including a dynamically consistent slip-rate function, 
applied to the 2019 Mw 7.1 Ridgecrest. Bull. Seism. Soc. Am., 112, 287-306. 
https://doi.org/10.1785/0120210138.

Skarlatoudis, A. A. et al. (2013). Ground-motion prediction equations 
of intermediate-depth earthquakes in the Hellenic arc, southern Aegean 
subduction area. Bull. Seism. Soc. Am., 103(3), 1952-1968. 
https://doi.org/10.1785/0120120265

------------------------------------------------------------------------
