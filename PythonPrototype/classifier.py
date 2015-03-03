# -*- coding: utf-8 -*-
"""
Created on Tue Mar 03 09:04:42 2015
Module for Feature extraction
@author: Patrick
"""

import math
import datastructure2 as ds

def testClassifier():
    data= ds.VTKData()
    data.loadFile('C:/out.1200.vtk')
    data.builtKDTree()
    
    classifier = Classifier(data)
    classifier.featureExtraction()
    
class Classifier:
    ##Class Variables
    _dataset = None
    _filePath = 'E:/modified.out.1200.vtk' 
    _featureExtractor = None
    _areaExtractor = None
    
    ##Class Methods
    def __init__(self, dataset):
        self._dataset = dataset
        self._featureExtractor = FeatureExtractor(self)
        self._areaExtractor = AreaExtractor(self)
        return
        
    def featureExtraction(self):
        """ Find features in whole dataset """
        ver_counter = 0
        for ver in self._dataset._vertexList:
            ver_counter +=1
            area = self._areaExtractor.getAreaNN(ver,6)
            featureValue = self._featureExtractor.VF_Angle_Diff(ver,area)
            ver._feature = featureValue
            if ((ver_counter*100.0/len(self._dataset._vertexList)) % 10)==0:
                        print('Feature computation reached: ' + str(ver_counter*100.0/len(self._dataset._vertexList)) + '%')
        self._dataset.writeFeatureToFile(self._filePath)
        return
        
class FeatureExtractor:
    ##Class Variables
    _classifier = None
    ##Class Methods
    def __init__(self,classifier):
        self._classifier = classifier
        return
    
    def VF_Angle_Diff(self,v0,area):
        """input v0: Vertex of Interest 
                area: tupple of list of distances and Verteces
           output : scalar value"""
        anglediff=[]        
        for ver in area[1]:
            #print(v0._ID,ver._ID)
            #print(ds.dot(v0._mag,ver._mag))
           # print(( v0._mag._length() * ver._mag._length() ) )
            diff =  math.acos(ds.dot(v0._mag,ver._mag) / ( v0._mag._length() * ver._mag._length() ) )
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
        
class AreaExtractor:
    ##Class Variables
    _classifier = None
    
    ##Class Methods
    def __init__(self,classifier):
        self._classifier = classifier
        return
    
    def getAreaNN(self,ver,n):
        """ input pos: Vertex of Interest  n: Number of nearest neighbours
            output: tupple list of distances and Verteces for specified Area definition"""
        xtupple = (ver._pos._x,ver._pos._y,ver._pos._z)
        d,ni = self._classifier._dataset._kdTree.query(xtupple,n)
  
        if ni == []:
            print("Cell not Found",ver._pos)
            return
        nn_list = []
        for i in ni:
            """remove duplicated IDS"""
            if i != ver._ID:
                nn_list.append(self._classifier._dataset._vertexList[i])   
        return d,nn_list
        
        
        
        
        