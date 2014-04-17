'''
This file is part of the lenstractor project.
Copyright 2012 David W. Hogg (NYU) and Phil Marshall (Oxford).

Description
-----------

Wrapper class for tractor operation, to enable high level options (ie sampling or optimization).

To-do
-----
- debug, test

'''

import numpy as np

import tractor
import lenstractor


# ============================================================================

class Explore():

      def __init__(self,by='sampling'):
            self.method = by

      def getName(self):
            return 'The Tractor`s Driver'

# ============================================================================

if __name__ == '__main__':

   pass
