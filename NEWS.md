# ESTER version 0.3.0

## Major updates

* ...
* removing distER.R and seqERboot.R (both deprecated in the last release)

## Minor updates

* doc updates
* remove verbose argument in simER.R
* adding the "correction" argument to aic.R for second-order bias correction

# ESTER version 0.2.0

## Major updates

* deprecating distER.R (will be removed in the next release)
* merging seqER.R and seqERboot.R into seqER.R, by adding a "nsims" argument in seqER.R
* deprecating seqERboot.R (will be removed in the next release)
* simER.R: parallelisation and update of the plot method
* introducing the analysER.R function which allow to analyse the results of simulations ran with simER.R

## Minor updates

* doc updates
* adding pseudo-BMA weights computations for brmsfit models in ictab.R
* seqER.R now accepts brmsfit models
* adding a "blind" argument to seqER.R to conduct triple blind analyses, with a "boundary" argument at which the sequential testing is stopped
