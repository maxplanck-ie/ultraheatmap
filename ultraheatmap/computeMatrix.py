import sys
import os
import numpy as np
import argparse

from deeptools.heatmapper import heatmapper
from deeptools.plotHeatmap import plotMatrix

sys.path.append(os.path.join(os.path.dirname(os.path.dirname((
                             os.path.realpath(__file__))))))


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
                  'proc number': configfile["numberOfProcessors"][index],
                  'sort regions': configfile["sortRegions"][index],
                  'sort using': configfile["sortUsing"][index],
                  'unscaled 5 prime': configfile["unscaled5prime"][index],
                  'unscaled 3 prime': configfile["unscaled3prime"][index]}
    return parameters


def __compute_matrix(regions, bigwigs, configfile, parameters, refIndex=None):
    """
    computing the corresponding matrix using deeptools/computeMatrix
    """
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
    hm.computeMatrix(score_file_list=bigwigs, regions_file=regions,
                     parameters=parameters,
                     blackListFileName=configfile["blackListFileName"][index],
                     verbose=configfile["verbose"][index], allArgs=matrix_args)

    return hm


def __clustering(hm, indexList, configfile):
    """
    """
    for index, i in enumerate(indexList):
        indexList[index] = i - 1
    if hm.parameters['min threshold'] is not None or\
            hm.parameters['max threshold'] is not None:
        hm.filterHeatmapValues(hm.parameters['min threshold'],
                               hm.parameters['max threshold'])

    if configfile["kmeans"] is not None:
        hm.matrix.hmcluster(configfile["kmeans"], method='kmeans')
    else:
        if configfile["hclust"] is not None:
            print("Performing hierarchical clustering."
                  "Please note that it might be very slow for large "
                  "datasets.\n")
            hm.matrix.hmcluster(configfile["hclust"], method='hierarchical')

    group_len_ratio = np.diff(hm.matrix.group_boundaries) /\
        len(hm.matrix.regions)
    if np.any(group_len_ratio < 5.0 / 1000):
        problem = np.flatnonzero(group_len_ratio < 5.0 / 1000)
        sys.stderr.write("WARNING: Group '{}' is too small for plotting,"
                         "you might want to remove it. "
                         "There will likely be an error message from "
                         "matplotlib regarding this "
                         "below.\n".format(hm.matrix.group_labels[problem[0]]))
    # TODO set sample & region labels!!
    if configfile["sortRegions"][0] != 'keep':
        hm.matrix.sort_groups(sort_using=configfile["sortUsing"][0],
                              sort_method=configfile["sortRegions"][0],
                              sample_list=indexList)
    outputMatrix_path = ""
    if configfile["outputReferenceMatrix"] is not None:
        outputMatrix_path = os.path.join(configfile["outputReferenceMatrix"])
        hm.save_matrix(outputMatrix_path)

    """ TODO: figure out how to do it directly from hm. For the moment when I
        use hm some parameters have incomapatible types. for example upstream
        is a value if it is read directly but is a list if it is read from a
        file."""
    if configfile["plotOutput"] is not None:
        if configfile["outputReferenceMatrix"] is None:
            outputMatrix_path = os.path.dirname(os.path.abspath(configfile["matrixOutput"]))
            outputMatrix_path += "/outputReferenceMatrix.gz"
            hm.save_matrix(outputMatrix_path)

        hm1 = heatmapper()
        hm1.read_matrix_file(outputMatrix_path)
        color_dict = {'colorMap': ['RdYlBu'],
                      'colorList': None,
                      'colorNumber': int(256),
                      'missingDataColor': 'black',
                      'alpha': float(1.0)}
        plotMatrix(hm1, os.path.join(configfile["plotOutput"]),
                   colorMapDict=color_dict)

    assert(configfile["outFileSortedRegions"])
    hm.save_BED(open(configfile["outFileSortedRegions"], "w"))


def sortbyreference(regions, bigwigs, indexList, configfile):
    # only on the ref.ones
    parameters = __parse_matrix_parameters(configfile, post_clustering=False)
    hm = __compute_matrix(regions, bigwigs, configfile, parameters,
                          refIndex=indexList)

    __clustering(hm, configfile["refIndices"], configfile)


def computefinalmatrix(regions, bigwigs, configfile):
    parameters = __parse_matrix_parameters(configfile, post_clustering=True)
    hm = __compute_matrix(regions, bigwigs, configfile, parameters,
                          refIndex=None)

    if configfile["samplesLabel"] and len(configfile["samplesLabel"]):
        hm.matrix.set_sample_labels(configfile["samplesLabel"])

    if configfile["regionsLabel"]:
        hm.matrix.set_group_labels(configfile["regionsLabel"])

    return hm
