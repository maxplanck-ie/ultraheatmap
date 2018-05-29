import pandas as pd
import csv
import os
import sys
from pybedtools import BedTool
import numpy as np
import gffutils
from deeptoolsapi.deeptoolsMatrix import Matrix


from coordinates.mapClosestGenes import keymap_from_closest_genes

def find_closest_genes(peaks, annotation, featureType, outputDir, filename = "test_mapped"):
    """
    Find the closest gene using bedtools.closest
    """
    Peaks = BedTool(peaks)
    Annotation = BedTool(annotation)

    Peaks=Peaks.sort()
    sites=Annotation.sort()
    if featureType:
        __filter_annotation(outputDir, featureType, annotation)
        sites=BedTool(outputDir+"filtered.gtf").sort()

    mapped=Peaks.closest(sites, t="first")

    if filename:
        mapped.saveas(filename)

    return(mapped)

def __filter_annotation(outputDir, featureType, annotation):
   """
   filters annotation for the gene type of the interest
   """
   with open(outputDir+"filtered.gtf","w") as filteredAnnotation:
        for feature in gffutils.DataIterator(annotation):
            if feature.featuretype == featureType:
               filteredAnnotation.write(str(feature)+'\n')

def extract_ge_folchange_per_peak(peaks, annotation, deseqtables, closestMapping,deseqfeature):
    """

    """
    #### <START>
    ## The following set of functions produce a keymap (dictionary) mapping
    ## each gene to associated peaks (<1:n>-mapping, n > 0)
    ##

    ## keymap: peak_key:gene_id
    ## gene_id
    ## peak_keys: cols 1-6 from bed format (1-based index)
    Peaks = BedTool(peaks)
    Peaks=Peaks.sort()
    keyMap_closest = keymap_from_closest_genes(closestMapping, Peaks)
    return(extractFoldChange(keyMap_closest, deseqtables,deseqfeature))

def __getValuesFromDEseqTable(geneid, deseqtable, deseqfeature):
    v = []
    for gid in geneid:
        if gid in deseqtable['GeneID'].values:
            x = float(deseqtable[deseqtable.GeneID == gid][deseqfeature])
            if np.isnan(x):
                x = np.nan
            v += [ x ]
        else:
            v += [ np.nan ]
    return v

def extractFoldChange(keyMap_closest, deseqtables, deseqfeature):
    """

    """
    matrixDict = {}
    geneIdtables =[parseGeneIdTable(table) for table in deseqtables]

    regions = [ key.split(';') for key in keyMap_closest ]

    valuesTab = np.empty((len(regions), len(deseqtables)), dtype=float)
    for i, table in enumerate(geneIdtables):
        values = __getValuesFromDEseqTable([keyMap_closest[key] for key in keyMap_closest], table, deseqfeature)
        valuesTab[:,i] = values


    return (Matrix(regions = regions, matrix = valuesTab, group_boundaries = [0,len(regions)], \
    sample_boundaries =  [x for x in range(0, len(deseqtables) + 1, 1)]))



def parseGeneIdTable(table_file):
    return(pd.read_csv(table_file, sep ='\t'))
