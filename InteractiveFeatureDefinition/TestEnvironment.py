# -*- coding: utf-8 -*-
"""
Created on Tue May 12 14:13:11 2015

@author: Patrick
"""
from ParaviewClassifier import *
from FeatureExtractorInterface import *
from FeatureExtractor_AngularDirectionChangingRate import *
from AreaExtractorSimple import *
from AreaExtractorInterface import *
from ClassifierInterface import *


def testClassifier(dataObject):
    classy = PWDatastructure(dataObject)
    areaEx = AreaExtractorSimple(classy)
    classy._areaExtractor = areaEx
    featureEx = FeatureExtratorADCR(classy,10,4) # check if magnetic field is number 8
    classy._featureExtractor = featureEx
    featureEx.evaluate()
    classy.addAttribute("adcr",1,featureEx._featureMap)
    

def main():
    """only dummy values here """


if __name__ == "__main__":
    main()