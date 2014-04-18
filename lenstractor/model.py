'''
This file is part of the lenstractor project.
Copyright 2012 David W. Hogg (NYU) and Phil Marshall (Oxford).

Description
-----------

Model class, to package up srcs and enable initialization.

To-do
-----
- debug, test

'''

import numpy as np

import tractor
import lenstractor


# ============================================================================

class Model():
    '''
    PURPOSE
      Define a named model that can be initialized, passed to a tractor etc.

    COMMENTS

    CONTENTS
      srcs             A list of tractor.sources
      name             etc

    OUTPUTS

    BUGS

    HISTORY
      2014-04-17       Started Marshall & Agnello (UCSB)
    '''
# ----------------------------------------------------------------------------
    
    def __init__(self,name,srcs=None):
    
        self.name = name
        
        if srcs == None:
            self.srcs = []
            self.initialize()
        
        return None
# ----------------------------------------------------------------------------
    
    def __str__(self):
        return '%s' % (self.name)

# ----------------------------------------------------------------------------
    
    def initialize(self,template=None):
                
        if template == None:
            pass
            
        else:
            
            pass
        
        return None
            
# ============================================================================

if __name__ == '__main__':

   pass
