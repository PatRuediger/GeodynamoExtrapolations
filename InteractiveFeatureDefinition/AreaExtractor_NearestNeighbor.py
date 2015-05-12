# -*- coding: utf-8 -*-
"""
Created on Tue May 12 13:18:31 2015

@author: Patrick

Area Extractor Implementation for Nearest Neighbor Area
"""

import numpy as np
from math import *
from abc import ABCMeta


class AreaExtractorNN(object):
    
    ##Class Members
    _classifier = None
    _VertexList = None  #stores all dataset vertex indices that are part of the area
    _DistanceList = None 
    
    ##Class Methods
    def __init__(self,classifier,center_pos,numNN):
        self._classifier = classifier
        self._DistanceList, self._VertexList = self._classifier.GetKDTree().query(center_pos,numNN) 
        
        
        
    
    
    def isPartofArea(self,pos):
        """ determines wheter a Point is part of the Area or not,
            espacially for points between two vertices in the _VertexList
            @input: numpy array
            @output: boolean
        """