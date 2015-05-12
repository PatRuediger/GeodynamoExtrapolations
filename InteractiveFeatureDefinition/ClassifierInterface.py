# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 10:17:03 2015

@author: Patrick Ruediger   
ruediger@rhrk.uni-kl.de

Interface for the use of the interactive feature defintion and extraction
"""

import numpy as np
from math import *
from abc import ABCMeta
from ParaviewClassifier import *


class ClassifierInterface:
    __metaclass__ = ABCMeta
    ##Class Members
    _dataSet =None
    _featureExtractor = None
    _areaExtractor = None

    #class methods
    def __init__(self,dataSet):
        self._dataSet = dataSet

    
    def bindExtrators(self):
        """ create and initialize the implemented Extrators
            and bind them to the appropriate variabels
        """
    
    
    
    def GetPoint(self,index):
        """ get Vertex Position by index
            @input: int
            @output: numpay array
        """
        
   
    def GetPointByNN(self,pos):
        """ get nearest Vertex Position in Dataset by a given input position
            @input: numpy array
            @output: numpy array
        """
        
  
    def GetPointDataByID(self,VertexID,ValueID):
        """ get Value of a Vertex by his ID and the ID of the Value if it has more then one
            -> support for multifield data
            @input: int, int
            @output: numpy array
        """
    
 
    def GetPointDataByPos(self,pos,ValueID):
        """ get Value of the nearest Vertex by the Position provided
            and the ID of the Value if it has more then one
            -> support for multifield data
            Which range metric is used, is defined in the implementation.
            Ensures the use of non-interpolated data points.
            @input: numpy array, int
            @output: numpy array
        """
         
    

    def addAttribute(self,name,numComp,indexArray,attArray):
        """ name : attrbute name
            numComp: 1 for sclars, 3 for vectors etc.
            indexArray: determins for which points you want to set the attribute
            attArray: list of attribute tuples to set, indexing should match the one from the indexArray
        """

  