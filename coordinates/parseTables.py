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
        assert filteredoutput
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

def extract_ge_folchange_per_peak(peaks, tables, closestMapping,features, IdColumn,hm):
    """

    """
    ## keyMap_closest: peak_key:gene_id
    ## peak_keys: cols 1-7 from bed format (1-based index)
    Peaks = BedTool(peaks)
    Peaks=Peaks.sort()
    keyMap_closest = keymap_from_closest_genes(closestMapping, Peaks)
    __update_matrix_values(peaks, keyMap_closest, tables,features,IdColumn,hm)
    __update_parameters(hm,len(tables))

def __getValuesFromGETable(peaks, keyMap_closest, table, features, IdColumn):
    print(features)
    v = []
    for peak in peaks:
        key = ';'.join(map(str,peak))
        value = keyMap_closest[key]
        if value in table[IdColumn].values:
            for feature in features: #TODO Check if pd.dataframe has an inbuilt function to get discontiniuos columns
                x = float(table[table[IdColumn] == value][feature])
                if np.isnan(x):
                     x = np.nan
                v += [ x ]*len(features)
        else:
            v += [ np.nan ]
    return v


def __getValuesFromTable(peaks, table, features, IdColumn):
    v = []
    for peak in peaks:
        name = peak[3] #TODO
        if name in table[IdColumn].values:
           for feature in features:
                x = float(table[table[IdColumn] == name][feature])
                if np.isnan(x):
                     x = np.nan
                v += [ x ]*len(features)
        else:
           v += [ np.nan ]
    return v

def __update_matrix_values(peaks, keyMap_closest, tables, features, IdColumn, hm):
    """

    """
    assert len(keyMap_closest) == len(peaks)
    valuesTab = np.empty((len(peaks), len(tables)), dtype=float)
    print(tables)
    for i, table in enumerate(tables):
        table = parseTable(table)
        values = __getValuesFromGETable(peaks, keyMap_closest, table, features, IdColumn)
        valuesTab[:,i] = values
        hm.matrix.sample_labels = hm.matrix.sample_labels + ["table"+str(i+1)]
    hm.matrix.matrix = np.concatenate((hm.matrix.matrix, valuesTab[:,]), axis = 1)
    last_col = hm.matrix.sample_boundaries[-1]
    hm.matrix.sample_boundaries = hm.matrix.sample_boundaries +[x+1+last_col for x in range(len(tables))]



def parseTable(table_file):
    return(pd.read_csv(table_file, sep ='\t'))

def parseMatrixRegions(regions):
    all_regions = []
    for group in  regions:
        for region in group:
            for start , end in region[1]:
                all_regions.append([region[0],start, end, region[2], region[3], region[4], region[5]])
    return all_regions

def __update_parameters(hm,length): ##XXX How????
    """

    """
    for i in range(length):
        hm.parameters['unscaled 5 prime'].append(0)
        hm.parameters['unscaled 3 prime'].append(0)
        hm.parameters['body'].append(1000)
        hm.parameters['downstream'].append(0)
        hm.parameters['upstream'].append(0)
        hm.parameters['ref point'].append(None)
        hm.parameters['bin size'].append(10)

def update_matrix_values(peaks, tables,features, IdColumn,hm):
    """
    The function is used for the one to one mapping between the name of the enriched regions
    from the given matrix and the name of the enriched regions from the provided tables. The
    correspoding values of each region, obtained from the tables, are added to the matrix values.
    """
    valuesTab = np.empty((len(peaks), len(tables)*len(features)), dtype=float)
    for i, table in enumerate(tables):
        table = parseTable(table)
        values = __getValuesFromTable(peaks, table, features, IdColumn)
        valuesTab[:,i] = values
        hm.matrix.sample_labels = hm.matrix.sample_labels + ["table"+str(i+1)]
    hm.matrix.matrix = np.concatenate((hm.matrix.matrix, valuesTab[:,]), axis = 1)
    last_col = hm.matrix.sample_boundaries[-1]
    hm.matrix.sample_boundaries = hm.matrix.sample_boundaries +[x+1+last_col for x in range(len(tables))]
    __update_parameters(hm,len(tables)*len(features))

def __read_tables_columns(tables, features):
    for table in tables:
         df = parseTable(table)
         for feature in features:
             if feature not in df.columns:
                sys.stderr.write("feature "+feature+" doesn't exist in table" + table) 
                exit(1)
         
   

