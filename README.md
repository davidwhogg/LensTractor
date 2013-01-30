## LensTractor

Finding gravitationally-lensed quasars in multi-epoch ground-based
imaging data, by model comparison with a flexible non-trivial galaxy
morphological model.

### Authors:

* Philip J. Marshall, University of Oxford
  <http://www.slac.stanford.edu/~pjm/>
* David W. Hogg, New York University
  <http://cosmo.nyu.edu/hogg/>
* Dustin Lang, Carnegie Mellon University

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
* `emcee`: As an alternative to optimization, LensTractor is MCMC-enabled, using Foreman-Mackey et al's emcee ensemble sampler.
* `astrometry.net`: LensTractor uses the `util` library to get its images' WCS right.

### Getting started:

Assuming you have numpy etc installed, first you'll need to download and install astrometry.net and The Tractor. 
Fortunately these come bundled together, so you can just do this:

 svn cat http://astrometry.net/svn/trunk/projects/tractor/checkout.sh | bash

Now check that The Tractor runs:

 cd tractor
 python tractor-sdss-synth.py

With the tractor on your python path you should be able to run the LensTractor example, as described in the output of 

 python LensTractor.py -h
