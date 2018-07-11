import pandas as pd
import csv
import os
import sys
from pybedtools import BedTool
import numpy as np
import gffutils

from coordinates.mapClosestGenes import keymap_from_closest_genes

def find_closest_genes(peaks, annotation, featureType, filteredoutput, filename = None):
    """
    Find the closest gene using bedtools.closest
    """
    Peaks = BedTool(peaks)
    Annotation = BedTool(annotation)
    Peaks=Peaks.sort()
    sites=Annotation.sort()
    if featureType:
        print(featureType)
        __filter_annotation(filteredoutput, featureType, annotation)
        sites=BedTool(filteredoutput).sort()
    mapped=Peaks.closest(sites, t="first")

    if filename:
        mapped.saveas(filename)

    return(mapped)

def __filter_annotation(filteredoutput, featureType, annotation):
   """
   filters annotation for the gene type of the interest
   """
   with open(filteredoutput,"w") as filteredAnnotation:
        for feature in gffutils.DataIterator(annotation):
            if feature.featuretype == featureType:
               filteredAnnotation.write(str(feature)+'\n')

def extract_ge_folchange_per_peak(peaks, tables, closestMapping,feature,IdColumn,hm):
    """

    """
    ## keyMap_closest: peak_key:gene_id
    ## peak_keys: cols 1-7 from bed format (1-based index)
    Peaks = BedTool(peaks)
    Peaks=Peaks.sort()
    keyMap_closest = keymap_from_closest_genes(closestMapping, Peaks)
    __extractFoldChange(peaks, keyMap_closest, tables,feature,IdColumn,hm)

def __getValuesFromTable(peaks, keyMap_closest, table, feature, IdColumn):
    v = []
    for peak in peaks:
        key = ';'.join(map(str,peak))
        print(key)
        value = keyMap_closest[key]
        if value in table[IdColumn].values:
            x = float(table[table[IdColumn] == value][feature])
            if np.isnan(x):
                x = np.nan
            v += [ x ]
        else:
            v += [ np.nan ]
    return v


def __extractFoldChange(peaks, keyMap_closest, tables, feature, IdColumn,hm):
    """

    """
    print(keyMap_closest)
    assert len(keyMap_closest) == len(peaks)
    valuesTab = np.empty((len(peaks), len(tables)), dtype=float)
    for i, table in enumerate(tables):
        table = parseTable(table)
        values = __getValuesFromTable(peaks, keyMap_closest, table, feature, IdColumn)
        valuesTab[:,i] = values
    hm.matrix.matrix = hm.matrix.matrix + valuesTab
    hm.matrix.sample_boundaries = hm.matrix.sample_boundaries +[x for x in range(0, len(tables) + 1, 1)]



def parseTable(table_file):
    return(pd.read_csv(table_file, sep ='\t'))

def parseMatrixRegions(regions):
    all_regions = []
    for group in  regions:
        for region in group:
            for start , end in region[1]:
                all_regions.append([region[0],start, end, region[2], region[3], region[4], region[5]])
    return all_regions
