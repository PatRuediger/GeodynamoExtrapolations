# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 12:16:42 2015

@author: Patrick Ruediger   
ruediger@rhrk.uni-kl.de

Interface for the use of the interactive feature defintion and extraction
"""

import numpy as np
from math import *
from abc import ABCMeta
from FeatureExtractor_AngularDirectionChangingRate import *


class FeatureExtratorInterface:
    __metaclass__ = ABCMeta
    
    ##Class Members
    _classifier = None
    """stores the evaluated features in a numpy array 
        with their respective VertexID
    """
    
    _featureMap = None 
    
     ##Class Methods
    def __init__(self,classifier):
        self._classifier = classifier
        
    
    def evaluate(self):
        """ evaluates the defined feature and fills the featureMap
        """
    

    def getInterpolatedFeature(self,pos):
        """ In Case the global interpolation method is not sufficient
            one can define a special one here
            @input: numpy array
            @output: numpy array
        """


        