## LensTractor

Finding gravitationally-lensed quasars in multi-epoch ground-based
imaging data, by model comparison with a flexible non-trivial galaxy
morphological model.

### Authors:

* Philip J. Marshall, Oxford University
  <http://www.slac.stanford.edu/~pjm/>
* David W. Hogg, New York University
  <http://cosmo.nyu.edu/hogg/>

### License:

Copyright 2012 Phil Marshall (Oxford) & David W. Hogg (NYU).  All
rights reserved.

### Dependencies:

* `numpy`
* `matplotlib`
* `Tractor`: LensTractor uses The Tractor (Hogg & Lang) to do the
  comparison between predicted image positions and the data: that is,
  The Tractor synthesizes images and computes likelihoods, and then 
  optimizes the model parameters.
