# ESTER version 0.2.0

## Major updates

* deprecating distER.R (will be removed in the next release)
* merging seqER and seqERboot into seqER, by adding a "nsims" argument in seqER
* deprecating seqERboot (will be removed in the next release)
* simER: parallelisation and update of the plot method, to plot distributions of ERs at each boundary (lower, upper, and nmax)
* introducing the compER function to compare the behaviour of different evidence ratios
* introducing the analysER function which allow to analyse the results of simulations ran with simER or compER

## Minor updates

* doc updates
* adding pseudo-BMA weights computations for brmsfit models in ictab.R
* adding a "blind" argument to seqER to conduct triple blind analyses, with a "bondary argument" at which the sequential testing is stopped

NB: we recommend using BIC-based rather than AIC-based evidence ratios for sequential testing
