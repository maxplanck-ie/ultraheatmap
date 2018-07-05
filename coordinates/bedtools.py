import pandas as pd
import csv
import os
import sys
from pybedtools import BedTool
import numpy as np
import gffutils
from deeptoolsapi.deeptoolsMatrix import Matrix


from coordinates.mapClosestGenes import keymap_from_closest_genes

def find_closest_genes(peaks, annotation, featureType, outputDir, filename = None):
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

def extract_ge_folchange_per_peak(peaks, deseqtables, closestMapping,deseqfeature,geneIdColumn):
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
    return(extractFoldChange(keyMap_closest, deseqtables,deseqfeature,geneIdColumn))

def __getValuesFromDEseqTable(geneid, deseqtable, deseqfeature, geneIdColumn):
    v = []
    for gid in geneid:
        if gid in deseqtable[geneIdColumn].values:
            x = float(deseqtable[deseqtable[geneIdColumn] == gid][deseqfeature])
            if np.isnan(x):
                x = np.nan
            v += [ x ]
        else:
            v += [ np.nan ]
    return v

def __parseRegions(keyMap_closest):
    """

    """
    regions =[]
    for key in keyMap_closest:
        region = key.split(';')
        chrom, start, end, name, score, strand = region[0:6]
        starts = start.split(",") ##XXX is it more than one?
        ends = end.split(",")
        regs = [(int(x), int(y)) for x, y in zip(starts, ends)]
        regions.append([chrom, regs, name, len(keyMap_closest), strand, score]) #XXX max_group_bound? I have set it the number of line since i am thinking that we always have one bed file at the time. Am I right?
    return regions


def extractFoldChange(keyMap_closest, deseqtables, deseqfeature, geneIdColumn):
    """

    """
    matrixDict = {}
    geneIdtables =[parseGeneIdTable(table) for table in deseqtables]

    valuesTab = np.empty((len(keyMap_closest), len(deseqtables)), dtype=float)
    for i, table in enumerate(geneIdtables):
        values = __getValuesFromDEseqTable([keyMap_closest[key] for key in keyMap_closest], table, deseqfeature, geneIdColumn)
        valuesTab[:,i] = values


    return (Matrix(regions = __parseRegions(keyMap_closest), matrix = valuesTab, group_boundaries = [0,len(keyMap_closest)], \
    sample_boundaries =  [x for x in range(0, len(deseqtables) + 1, 1)]))



def parseGeneIdTable(table_file):
    return(pd.read_csv(table_file, sep ='\t'))
