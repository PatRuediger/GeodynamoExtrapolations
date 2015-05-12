# -*- coding: utf-8 -*-
"""
Created on Tue May 12 14:29:10 2015

@author: Patrick
"""

import numpy as np
from math import *
from abc import ABCMeta


class AreaExtractorSimple(object):
    
    ##Class Members
    _classifier = None
    _VertexList = None  #stores all dataset vertex indices that are part of the area
    _DistanceList = None 
    
    ##Class Methods
    def __init__(self,classifier):
        self._classifier = classifier
        self._VertexList = self._classifier._dataRef.GetPoints() 
        
        
        
    
    
    def isPartofArea(self,pos):
        """ determines wheter a Point is part of the Area or not,
            espacially for points between two vertices in the _VertexList
            @input: numpy array
            @output: boolean
        """