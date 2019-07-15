#!/usr/bin/env python

import pandas as pd
import csv
import os
import sys
from pybedtools import BedTool
import numpy as np
import gffutils

from ultraheatmap.mapClosestGenes import keymap_from_closest_genes

# TODO:
def find_closest_genes(peaks, annotation, annotationFeature, filteredoutput,
                       referencePoint, filename=None):
    """
    Find the closest gene using bedtools.closest
    """
    Peaks = BedTool(peaks)
    Annotation = BedTool(annotation)
    Peaks = Peaks.sort()
    sites = Annotation.sort()
    if annotationFeature:
        filteredAnnotation = __filter_annotation(filteredoutput,
                                                 annotationFeature, annotation,
                                                 referencePoint)
        sites = BedTool(filteredAnnotation).sort()
    elif referencePoint:
        filteredAnnotation = list()
        for feature in gffutils.DataIterator(annotation):
            filteredAnnotation.append(str(__get_reference_coordinate(feature,
                                          referencePoint)))
        sites = BedTool(filteredAnnotation).sort()
    mapped = Peaks.closest(sites, t="first")

    if filename:
        mapped.saveas(filename)

    return(mapped)


def __get_reference_coordinate(feature, referencePoint):
    """
    """
    if (feature.strand == '+' and referencePoint == 'TSS') or\
       (feature.strand == '-' and referencePoint == 'TES'):
        feature.end = feature.start
    else:
        assert(feature.strand == '-' and referencePoint == 'TSS') or\
              (feature.strand == '+' and referencePoint == 'TES')
        feature.start = feature.end
    return feature


# TODO
def __filter_annotation(filteredoutput, annotationFeature, annotation,
                        referencePoint):
    """
    Filter the annotation for the gene type of interest.
    """
    if filteredoutput:
        fo = open(filteredoutput, "w")
    filteredAnnotation = list()
    for feature in gffutils.DataIterator(annotation):
        if feature.featuretype == annotationFeature:
            if filteredoutput:
                fo.write(str(feature)+'\n')
            if referencePoint:
                feature = __get_reference_coordinate(feature, referencePoint)
            filteredAnnotation.append(str(feature))

    return filteredAnnotation


def extract_ge_folchange_per_peak(peaks, tables, closestMapping, features,
                                  IdColumn, hm):
    """

    """
    ## keyMap_closest: peak_key:gene_id
    ## peak_keys: cols 1-7 from bed format (1-based index)
    Peaks = BedTool(peaks)
    Peaks=Peaks.sort()
    keyMap_closest = keymap_from_closest_genes(closestMapping, Peaks)
    __update_matrix_values(peaks, keyMap_closest, tables,features,IdColumn,hm)

def __getValuesFromGETable(peaks, keyMap_closest, table, features, IdColumn):
    v = np.empty((len(peaks), len(features)), dtype=float)
    for i, peak in enumerate(peaks):
        key = ';'.join(map(str,peak))
        value = keyMap_closest[key]
        if value in table[IdColumn].values: #value is geneId
            for j, feature in enumerate(features): #TODO Check if pd.dataframe has an inbuilt function to get discontiniuos columns
                x = float(table[table[IdColumn] == value][feature])
                if np.isnan(x):
                     x = np.nan
                v[i,j] = x
        else:
            v[i] = [ np.nan ]*len(features)
    return v


def __getValuesFromNameTable(peaks, table, features, IdColumn):
    v = np.empty((len(peaks), len(features)), dtype=float)
    for i, peak in enumerate(peaks):
        name = peak[3]
        if name in table[IdColumn].values:
           for j, feature in enumerate(features):
                x = float(table[table[IdColumn] == name][feature])
                if np.isnan(x):
                     x = np.nan
                v[i,j] = x
        else:
           v[i] = [ np.nan ]*len(features)
    return v

def __update_matrix_values(peaks, keyMap_closest, tables, features, IdColumn, hm): #TODO two different update_value function can be merged into one, just need to an if / else over the arg.annotation
    """

    """
    assert len(keyMap_closest) == len(peaks)
    valuesTab = np.empty((len(peaks), len(tables)*len(features)), dtype=float)
    print(tables)
    for i, table in enumerate(tables):
        table = parseTable(table)
        values = __getValuesFromGETable(peaks, keyMap_closest, table, features, IdColumn)
        valuesTab[:,i*len(features):(i*len(features)+len(features))] = values
        for feature in features:
            hm.matrix.sample_labels = hm.matrix.sample_labels + ["table"+str(i)+"_"+feature]
    hm.matrix.matrix = np.concatenate((hm.matrix.matrix, valuesTab[:,]), axis = 1)
    current_last_col = hm.matrix.sample_boundaries[-1]
    hm.matrix.sample_boundaries = hm.matrix.sample_boundaries +[x+1+current_last_col for x in range(len(tables)*len(features))]
    __update_parameters(hm,len(tables)*len(features))




def parseTable(table_file):
    return(pd.read_csv(table_file, sep ='\t'))

def parseMatrixRegions(regions):
    all_regions = []
    for group in  regions:
        for region in group:
            for start , end in region[1]:
                all_regions.append([region[0],start, end, region[2], region[3], region[4], region[5]])
    return all_regions

def __update_parameters(hm,length):
    """

    """
    print(length)
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
    print(valuesTab.shape)
    for i, table in enumerate(tables):
        table = parseTable(table)
        values = __getValuesFromNameTable(peaks, table, features, IdColumn)
        valuesTab[:,i*len(features):(i*len(features)+len(features))] = values
        for feature in features:
            hm.matrix.sample_labels = hm.matrix.sample_labels + ["table"+str(i)+"_"+feature]
    hm.matrix.matrix = np.concatenate((hm.matrix.matrix, valuesTab[:,]), axis = 1)
    current_last_col = hm.matrix.sample_boundaries[-1]
    hm.matrix.sample_boundaries = hm.matrix.sample_boundaries +[x+1+current_last_col for x in range(len(tables)*len(features))]
    __update_parameters(hm,len(tables)*len(features))

def __read_tables_columns(tables, features):
    for table in tables:
         df = parseTable(table)
         for feature in features:
             if feature not in df.columns:
                sys.stderr.write("feature "+feature+" doesn't exist in table" + table)
                exit(1)
