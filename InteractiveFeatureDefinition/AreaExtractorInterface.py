# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 12:01:26 2015

@author: Patrick Ruediger   
ruediger@rhrk.uni-kl.de

Interface for the use of the interactive feature defintion and extraction
"""

import numpy as np
from math import *
from abc import ABCMeta
from AreaExtractor_NearestNeighbor import *
from AreaExtractorSimple import *


class AreaExtractorInterface:
    __metaclass__ = ABCMeta
    
    ##Class Members
    _classifier = None
    _VertexList = None  #stores all dataset vertex indices that are part of the area
    
    ##Class Methods
    def __init__(self,classifier):
        self._classifier = classifier
    
    def isPartofArea(self,pos):
        """ determines wheter a Point is part of the Area or not,
            espacially for points between two vertices in the _VertexList
            @input: numpy array
            @output: boolean
        """
    
