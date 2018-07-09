import sys
import os
import pandas as pd
import numpy as np

sys.path.append(os.path.join(os.path.dirname(os.path.dirname((os.path.realpath(__file__))))))
from deeptoolsapi.plotHeatMap import __plot_heatmap
from deeptools.heatmapper import heatmapper

def __compute_matrix(bw, bed, configfile): ##XXX What about metagene?
   """
   computing the corresponding matrix using deeptools/computeMatrix
   """
   parameters = {'upstream': configfile["beforeRegionStartLength"],
                  'downstream': configfile["afterRegionStartLength"],
                  'body': configfile["regionBodyLength"],
                  'bin size': configfile["binSize"],
                  'ref point': configfile["referencePoint"],
                  'verbose': configfile["verbose"],
                  'bin avg type': configfile["averageTypeBins"],
                  'missing data as zero': configfile["missingDataAsZero"],
                  'min threshold': configfile["minThreshold"],
                  'max threshold': configfile["maxThreshold"],
                  'scale': configfile["scale"],
                  'skip zeros': configfile["skipZeros"],
                  'nan after end': configfile["nanAfterEnd"],
                  'proc number': configfile["numberOfProcessors"],
                  'sort regions': configfile["sortRegions"],
                  'sort using': configfile["sortUsing"],
                  'unscaled 5 prime': configfile["unscaled5prime"],
                  'unscaled 3 prime': configfile["unscaled3prime"]
   }
   hm = heatmapper()
   hm.computeMatrix(score_file_list = bw, regions_file = bed, parameters = parameters, blackListFileName=None, verbose=False, allArgs=None)
   return hm

def sortbyreference(regions,refIndex,bigwigs,configfile):
    refList=[]
    for index in refIndex:
        refList.append(bigwigs[int(index)])
    #only on the ref.ones
    ref_bw=" ".join(refList)
    hm = __compute_matrix(refList, regions, configfile)
    __plot_heatmap(hm, refIndex, configfile)


def computefinalmatrix(regions, bigwigs, configfile):
    hm = __compute_matrix(bigwigs, regions, configfile)
    return hm
