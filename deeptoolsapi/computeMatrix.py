import sys
import os
import pandas as pd
import numpy as np

sys.path.append(os.path.join(os.path.dirname(os.path.dirname((os.path.realpath(__file__))))))
from deeptoolsapi.plotHeatMap import __plot_heatmap
from deeptools.heatmapper import heatmapper

def __compute_matrix(bw, bed, configfile, args, pre_cluster_mode, boundries): ##XXX What about metagene?
   """
   computing the corresponding matrix using deeptools/computeMatrix
   """
   region_body =""
   upstream = None
   downstream = None
   if pre_cluster_mode == 'reference-point':
       region_body = 0
       assert(len(boundries)==2)
       upstream = boundries[0]
       downstream = boundries[1]
   elif pre_cluster_mode == 'scale-region': #TODO make them user defined values
       region_body = 1000
       upstream = 0
       downstream = 0
   else:
       print(pre_cluster_mode)
       assert(pre_cluster_mode=="")
       assert(boundries==[])
       region_body = configfile["regionBodyLength"]
       upstream = configfile["beforeRegionStartLength"]
       downstream = configfile["afterRegionStartLength"]
   parameters = {'upstream': upstream,
                  'downstream': downstream,
                  'body': region_body,
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
   hm.computeMatrix(score_file_list = bw, regions_file = bed, parameters = parameters, blackListFileName=None, verbose=False, allArgs=args)
   return hm

def sortbyreference(regions,refIndex,bigwigs,configfile, args, pre_cluster_mode, boundries):
    refList=[]
    for index in refIndex:
        assert(int(index) >= 1)
        refList.append(bigwigs[int(index)-1])
    #only on the ref.ones
    ref_bw=" ".join(refList)
    hm = __compute_matrix(refList, regions, configfile, args, pre_cluster_mode, boundries)
    __plot_heatmap(hm, refIndex, configfile)


def computefinalmatrix(regions, bigwigs, configfile, args):
    pre_cluster_mode =""
    boundries=[]
    hm = __compute_matrix(bigwigs, regions, configfile, args, pre_cluster_mode, boundries)
    return hm
