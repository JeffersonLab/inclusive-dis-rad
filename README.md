# inclusive-dis-rad
inclusive DIS generator, taking as input Bosted and CJ15 parametrizations with radiative corrections, can generate and work with grids 

To install `git clone https://github.com/JeffersonLab/inclusive-dis-rad.git`

cd inclusive-dis-rad make

Need to define the path for the pdf sets using env variable DISRAD_PDF to define the path to SF grid, when running with --xgrid option

Example `./generate-dis --trig 10000000 --nmax 10000` will write 100 gemc lund type data files with 10K DIS electrons in the current directory



Example 2 `./generate-dis --trig 10000 --docker` will write a single gemc lund type data file dis-rad.dat with 10K events  

Example 3 `./generate-dis --trig 10000 --rad 3` will produce also a radiative photon and store the RC factor in the header of the event

To get more command line options
`./generate-dis --help`
