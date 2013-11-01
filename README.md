## LensTractor

Finding gravitationally-lensed quasars in multi-epoch ground-based
imaging data, by model comparison with a flexible non-trivial galaxy
morphological model.

### Authors:

* Philip J. Marshall, KIPAC, SLAC
  <http://www.slac.stanford.edu/~pjm/>
* David W. Hogg, CCPP, New York University
  <http://cosmo.nyu.edu/hogg/>
* Dustin Lang, Physics dept, Carnegie Mellon University


### License:

LensTractor is being developed by Phil Marshall (KIPAC), David W. Hogg (NYU) and Dustin Lang (CMU). We welcome collaborators: if you are interested in using and improving LensTractor, please get in touch with Phil at `pjm at slac.stanford.edu`. The code is distributed under the GPL v2 license, which means that you can copy and edit it however you like, as long as you also distribute the result under the GPL v2 license. In practice, we hope that you will fork the code on github (from where you can easily ingest the updates we make, and from where we can pull in the excellent changes you make!) and keep in touch, so we can make sure that everyone's work is properly accredited in the publications that arise from LensTractor's testing and application. 

At the very least, we ask you to cite this website when either making use of the code, or the ideas and algorithms it contains.


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

With the tractor on your python path you should be able to run the LensTractor example, as described in LensTractor.py. For help, see the output of 

    python LensTractor.py -h
