# ============================================================================
'''
This file is part of the lenstractor project.
Copyright 2012 David W. Hogg (NYU) and Phil Marshall (Oxford).

Description
-----------

Plots made by lenstractor, including: 
* Progress on all images 

'''

import pylab as plt

import numpy as np

import tractor
import lenstractor

# Global variables:
px,py = 2,2
figprops = dict(figsize=(5*px,5*py), dpi=128)
adjustprops = dict(\
   left=0.05,\
   bottom=0.05,\
   right=0.95,\
   top=0.95,\
   wspace=0.1,\
   hspace=0.1)

# ============================================================================
# All progress plots.

def Plot_state(t,suffix):
   '''
   Make all the plots we need to assess the state of the tractor. Mainly, 
   a multi-panel figure of image, synthetic image and chi, for each image being 
   modelled.
   
   t is a Tractor object, containing a list of images.
   '''
   
   for i,image in enumerate(t.images):

      if image.name is None:
         imname = suffix+str(i)
      else:
         imname = image.name

      scale = np.sqrt(np.median(1.0/image.invvar[image.invvar > 0.0]))
 
      ima = dict(interpolation='nearest', origin='lower',
                     vmin=-30.*scale, vmax=3.*scale)
      chia = dict(interpolation='nearest', origin='lower',
                        vmin=-5., vmax=5.)
      psfa = dict(interpolation='nearest', origin='lower')

      fig = plt.figure(**figprops)
      fig.subplots_adjust(**adjustprops)
      plt.clf()
      plt.gray()

      plt.subplot(py,px,1)
      plt.imshow(-image.data, **ima)
      tidyup_plot()
      plt.title('Observed image')

      model = t.getModelImages()[i]
      # print "lenstractor.Plot_state: minmax of model = ",np.min(model),np.max(model)
      plt.subplot(py,px,2)
      plt.imshow(-model, **ima)
      tidyup_plot()
      plt.title('Predicted image')

      chi = t.getChiImage(i)
      plt.subplot(py,px,3)
      plt.imshow(-chi, **chia)
      tidyup_plot()
      plt.title('Residuals ($\pm 5\sigma$)')

      psfimage = image.psf.getPointSourcePatch(*model.shape).patch
      plt.subplot(py,px,4)
      plt.imshow(-psfimage, **psfa)
      tidyup_plot()
      plt.title('PSF')

      plt.savefig(imname+'_'+suffix+'.png')   
   
   return

# ----------------------------------------------------------------------------
# Turn off the axis ticks and labels:

def tidyup_plot():
   ax = plt.gca()
   ax.xaxis.set_ticks([])
   ax.yaxis.set_ticks([])
   return

# ============================================================================
