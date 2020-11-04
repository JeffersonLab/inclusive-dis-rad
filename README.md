# inclusive-dis-rad
inclusive DIS generator, taking as input Bosted and CJ15 parametrizations with radiative corrections, can generate and work with grids 

To install `git clone https://github.com/JeffersonLab/inclusive-dis-rad.git`

cd inclusive-dis-rad make

Need to define the path for the pdf sets using env variable DISRAD_PDF to define the path to SF grid, when running with --xgrid option

Example `./inclusive-dis-rad --trig 10000000 --nmax 10000` will write 100 gemc lund type data files with 10K DIS electrons in the current directory



Example 2 `./inclusive-dis-rad --trig 10000 --docker` will write a single gemc lund type data file dis-rad.dat with 10K events  

Example 3 `./inclusive-dis-rad --trig 10000 --rad 3` will produce also a radiative photon and store the RC factor in the header of the event

To get more command line options
`./inclusive-dis-rad --help`

To get a new random seed use in command line options
`./inclusive-dis-rad --seed NNNN`


Default option uses Bosted's parameterization coded in F1F209.f

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c F1F209.f
c Package of FORTRAN subroutine describing fits to inclusive inelastic
c electron scattering from proton, neutron, deuteron, or heaviere nuclei
c Proton fit is described in:
c   M.E. Christy and P.E. Bosted, ``Empirical Fit to Precision 
c    Inclusive Electron-Proton Cross Sections in the Resonance Region'',
c    (arXiv:0712.3731). Submitted to Phys. Rev. C.
c Deuteron (and netron) fit is described in:
c    P.E. Bosted and M.E. Christy, ``Empirical Fit to Inelastic 
c    Electron-Deuteron and Electron-Neutron
c    Resonance Region Transverse Cross Sections, 
c    (arXiv:0711.0159), publichished in Phys. Rev. C 77, 065206 (2008). (
c New fits for A>2 by Vahe M. and Peter B. (to be publsihed).
