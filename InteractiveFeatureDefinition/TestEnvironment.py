# -*- coding: utf-8 -*-
"""
Created on Tue May 12 14:13:11 2015

@author: Patrick
"""

import sys
sys.path.append('/Users/patrickrudiger/GeodynamoExtrapolations/InteractiveFeatureDefinition')

from ParaviewClassifier import *
from FeatureExtractorInterface import *
from FeatureExtractor_AngularDirectionChangingRate import *
from AreaExtractorSimple import *
from AreaExtractorInterface import *
from ClassifierInterface import *



def testClassifier(dataObject,dataOutput):
    classy = PWDatastructure(dataObject)
    areaEx = AreaExtractorSimple(classy)
    classy._areaExtractor = areaEx
    featureEx = FeatureExtratorADCR(classy,10,"magnetic_field") 
    classy._featureExtractor = featureEx
    featureEx.evaluate()
    classy.addAttribute("adcr",1,featureEx._featureMap,dataOutput)
    

def main():
	import sys
	sys.path.append('/Users/patrickrudiger/GeodynamoExtrapolations/InteractiveFeatureDefinition')
	import paraview
	import TestEnvironment as TE
	data = self.GetInput()
	dataOut = self.GetOutput()
	TE.testClassifier(data,dataOut)

if __name__ == "__main__":
    main()