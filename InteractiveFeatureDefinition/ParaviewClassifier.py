# -*- coding: utf-8 -*-
"""
Created on Fri May 08 10:57:03 2015

@author: Patrick Ruediger   
ruediger@rhrk.uni-kl.de

Interface for the use of the interactive feature defintion and extraction
"""
import numpy as np
from math import *
from abc import ABCMeta
from ClassifierInterface import *
#from paraview.simple import *
import vtk
from scipy.spatial import KDTree

class PWDatastructure(object):
    """ A simple datastructe which maps between the vtk data output format and 
        an easy to use format for feature extraction.
    """
    _dataSet =None
    _featureExtractor = None
    _areaExtractor = None
    _dataRef = None
    _numPts = None
    _dataRef_KDTree = None    
        
    def __init__(self, dataRef):
        self._dataRef = dataRef   
        self._numPts = self._dataRef.GetPoints().GetNumberOfPoints()
        print "Preocessing " + str(self._numPts) + " Data Points"
        vertexListTriples = []
        for i in range(0,self._numPts):
            vertexListTriples.append(self._dataRef.GetPoint(i))
        self._dataRef_KDTree = KDTree(vertexListTriples)

    #to be implemented from interface
    def GetPoint(self,index):
        """returns an array of doubles"""
        return self._dataRef.GetPoint(index)
        
    def GetPointDataByName(self,name,index):
        """ name: should be the exact name of the attribute as in the data set
            index:  the refering vertex/point index for which you want the attribute
            returns a Tupple if not a scalar (else a double)"""
        output = self._dataRef.GetPointData().GetArray(name).GetTuple3(index)
        #if (output[1]<1.0e-200) or (output[2]<1.0e-200):
        #    return output[1]
        #else : 
        return output
            
    #to be implemented from interface        
    def GetPointDataByID(self,ID,index):    
        """ name: should be the exact name of the attribute as in the data set
            index:  the refering vertex/point index for which you want the attribute
            returns a Tupple if not a scalar (else a double)"""
        output = self._dataRef.GetPointData().GetArray(ID).GetTuple3(index)
        #if (output[1]<1.0e-200) or (output[2]<1.0e-200):
        #    return output[1]
        #else : 
        return output
            
    #to be implemented from interface
    def GetPointByNN(self,pos):
        """ get nearest Vertex Position in Dataset by a given input position
            @input: numpy array
            @output: numpy array
        """
        d, nn = self._dataRef_KDTree.query(pos,1)
        return self.GetPoint(nn)            

     #to be implemented from interface
    def GetPointDataByPos(self,pos,ValueID):
        """ get Value of the nearest Vertex by the Position provided
            and the ID of the Value if it has more then one
            -> support for multifield data
            Which range metric is used, is defined in the implementation.
            Ensures the use of non-interpolated data points.
            @input: numpy array, int
            @output: numpy array
        """
        d, nn = self._dataRef_KDTree.query(pos,1)
        return self.GetPointDataByID(nn,ValueID)
        

     #to be implemented from interface
    def addAttribute(self,name,numComp,featureMap,dataOutput):
        """ name : attrbute name
            numComp: 1 for sclars, 3 for vectors etc.
            indexArray: determins for which points you want to set the attribute
            attArray: list of attribute tuples to set, indexing should match the one from the indexArray
            """
        ca = vtk.vtkFloatArray()
        ca.SetName(name)
        ca.SetNumberOfComponents(numComp) #(for scalars 1)
        ca.SetNumberOfTuples(self._numPts)
        
        for i in range(0,len(featureMap)):
            if numComp == 1:
                ca.SetTuple1(featureMap.keys()[i],featureMap.values()[i])
            if numComp == 3:
                ca.SetTuple3(featureMap.keys()[i],featureMap.values()[i])
        # if additional ones are needed just add the cases here, see vtkDataArray::setTuple for possible function calls
        dataOutput.GetPointData().AddArray(ca)
            
    def GetKDTree(self):
        return self._dataRef_KDTree