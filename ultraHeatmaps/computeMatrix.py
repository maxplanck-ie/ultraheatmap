import sys
import os
import pandas as pd
import numpy as np
import argparse

sys.path.append(os.path.join(os.path.dirname(os.path.dirname((os.path.realpath(__file__))))))
from deeptools.heatmapper import heatmapper

def __parse_matrix_parameters(configfile, post_clustering=False):
    index = 0
    if post_clustering:
        index = 1
    minThreshold = 0
    maxThreshold = 0
    if configfile["minThreshold"][index] == '':
        minThreshold = None
    else:
        minThreshold = configfile["minThreshold"][index]
    if configfile["maxThreshold"][index] == '':
        maxThreshold = None
    else:
        maxThreshold = configfile["maxThreshold"][index]
    parameters = {'upstream': configfile["beforeRegionStartLength"][index],
                   'downstream': configfile["afterRegionStartLength"][index],
                   'body': configfile["regionBodyLength"][index],
                   'bin size': configfile["binSize"][index],
                   'ref point': configfile["referencePoint"][index],
                   'verbose': configfile["verbose"][index],
                   'bin avg type': configfile["averageTypeBins"][index],
                   'missing data as zero': configfile["missingDataAsZero"][index],
                   'min threshold': minThreshold,
                   'max threshold': maxThreshold,
                   'scale': configfile["scale"][index],
                   'skip zeros': configfile["skipZeros"][index],
                   'nan after end': configfile["nanAfterEnd"][index],
                   'proc number': configfile["numberOfProcessors"], # some for both runs
                   'sort regions': configfile["sortRegions"][index],
                   'sort using': configfile["sortUsing"][index],
                   'unscaled 5 prime': configfile["unscaled5prime"][index],
                   'unscaled 3 prime': configfile["unscaled3prime"][index]
    }
    return parameters

def __compute_matrix(regions, bigwigs, configfile, parameters, refIndex = None):
   """
   computing the corresponding matrix using deeptools/computeMatrix
   """
<<<<<<< HEAD
   region_body =""
   upstream = None
   downstream = None
   referencePoint = None
   if pre_cluster_mode == 'reference-point':
       region_body = 0
       assert(len(boundries)==2)
       upstream = boundries[0]
       downstream = boundries[1]
       referencePoint = 'TSS'
   elif pre_cluster_mode == 'scale-regions':
       region_body = 1000
       upstream = 0
       downstream = 0
   else:
       assert(pre_cluster_mode=="")
       assert(boundries==[])
       region_body = configfile["regionBodyLength"]
       upstream = configfile["beforeRegionStartLength"]
       downstream = configfile["afterRegionStartLength"]
       referencePoint = configfile["referencePoint"]

   parameters = {'upstream': upstream,
                  'downstream': downstream,
                  'body': region_body,
                  'bin size': configfile["binSize"],
                  'ref point': referencePoint,
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
=======
>>>>>>> doc_new
   hm = heatmapper()

   if refIndex:
       bigwigs = [bigwigs[int(i)-1] for i in refIndex]

   index = 1 if not refIndex else 0

   matrix_args = argparse.Namespace()
   matrix_args.transcriptID = configfile['transcriptID'][index]
   matrix_args.exonID = configfile['exonID'][index]
   matrix_args.transcript_id_designator = configfile['transcript_id_designator'][index]
   matrix_args.samplesLabel = configfile['samplesLabel']
   matrix_args.exonID = configfile['exonID'][index]

   hm.computeMatrix(score_file_list = bigwigs, regions_file = regions, parameters = parameters,
    blackListFileName=configfile["blackListFileName"][index], verbose=configfile["verbose"][index], allArgs=matrix_args)

   return hm

<<<<<<< HEAD
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
=======

def __clustering(hm,indexList, configfile):
   """
   """
   if hm.parameters['min threshold'] is not None or hm.parameters['max threshold'] is not None:
        hm.filterHeatmapValues(hm.parameters['min threshold'], hm.parameters['max threshold'])

   if configfile["sortRegions"][0] == 'keep':
      configfile["sortRegions"][0] = 'no'  # These are the same thing ##XXX what does that mean???

   if configfile["kmeans"] is not None:
        hm.matrix.hmcluster(configfile["kmeans"], method='kmeans')
   else:
        if configfile["hclust"] is not None:
            print("Performing hierarchical clustering."
                  "Please note that it might be very slow for large datasets.\n")
            hm.matrix.hmcluster(configfile["hclust"], method='hierarchical')

   group_len_ratio = np.diff(hm.matrix.group_boundaries) / len(hm.matrix.regions)
   if np.any(group_len_ratio < 5.0 / 1000):
        problem = np.flatnonzero(group_len_ratio < 5.0 / 1000)
        sys.stderr.write("WARNING: Group '{}' is too small for plotting, you might want to remove it. "
                         "There will likely be an error message from matplotlib regarding this "
                         "below.\n".format(hm.matrix.group_labels[problem[0]]))

   if configfile["sortRegions"][0] != 'no':
        hm.matrix.sort_groups(sort_using=configfile["sortUsing"][0],sort_method=configfile["sortRegions"][0],sample_list=indexList)

   assert(configfile["outFileSortedRegions"])
   hm.save_BED(open(configfile["outFileSortedRegions"], "w"))


def sortbyreference(regions, bigwigs, indexList, configfile):
    # only on the ref.ones
    parameters = __parse_matrix_parameters(configfile, post_clustering = False)
    hm = __compute_matrix(regions, bigwigs, configfile, parameters, refIndex = indexList)
    __clustering(hm, configfile["refIndex"], configfile)


def computefinalmatrix(regions, bigwigs, configfile):
    parameters = __parse_matrix_parameters(configfile, post_clustering = True)
    hm = __compute_matrix(regions, bigwigs, configfile, parameters, refIndex = None)
>>>>>>> doc_new

    if configfile["samplesLabel"] and len(configfile["samplesLabel"]):
         hm.matrix.set_sample_labels(args.samplesLabel)

    if configfile["regionsLabel"]:
         hm.matrix.set_group_labels(configfile["regionsLabel"])


    return hm
