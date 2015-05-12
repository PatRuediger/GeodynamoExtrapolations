# -*- coding: utf-8 -*-
"""
Created on Tue May 12 13:34:31 2015

@author: Patrick
Implemented Feature Extractor for Angular Direction Changing Rate

"""

import numpy as np
from math import *
from abc import ABCMeta

class FeatureExtratorADCR(object):

    
    ##Class Members
    _classifier = None
    """stores the evaluated features in an array 
        with their respective VertexID
    """
    _NNRange = 0
    _featureMap = dict() 
    _valueID = 0
    
     ##Class Methods
    def __init__(self,classifier,NNRange,valueID):
        self._classifier = classifier
        self._NNRange = NNRange
        self._valueID =valueID

        

    def VF_Angle_Diff(self,v0,area):
        """input v0: Vertex ID of Interest 
                area: tupple of list of distances and Verteces
           output : scalar value"""
        anglediff=[]        
      #  v0p = self._classifier.GetPoint(v0)
        v0v = self._classifier.GetPointDataByID(v0,self._valueID)
        for ver in area[1]:
           # verp =  self._classifier.GetPoint(ver)
            verv =  self._classifier.GetPointDataByID(ver,self._valueID)
            diff =  math.acos(np.dot(v0v,verv) / (np.linalg.norm(v0v)*np.linalg.norm(verv) ) )
            anglediff.append(diff)
        output = self.getFeatureValue(anglediff,area[0])
        return output
        
    def getFeatureValue(self,diffs,weights):
        """input diffs: list of feature differences in area
                 weights: list of weightening facors
            output scalar value... may be changed later"""
        #unweighted for testing
        weightedValues = []
        for diff in diffs:
            if diff < math.pi :
                diffValue = math.exp(0.5*(diff- math.pi))
                weightedValues.append(diffValue)
            else:
                diffValue = math.exp(0.5*(-diff + math.pi))
                weightedValues.append(diffValue)
        mean = sum(weightedValues)/len(weightedValues)
        return mean      
        
    
    def evaluate(self):
        """ evaluates the defined feature and fills the featureMap
        """
        area = self._classifier._areaExtractor._VertexList , self._classifier._areaExtractor._DistanceList
        for ver in range(0,self._classifier._numPts):
            subArea = self._classifier.GetKDTree().query(area[0].GetPoint(ver),NNRange)
            value = self.VF_Angle_Diff(ver,subArea)
            self._featureMap[ver] = value
    
    def getInterpolatedFeature(self,pos):
        """ In Case the global interpolation method is not sufficient
            one can define a special one here
            @input: numpy array
            @output: numpy array
        """
        