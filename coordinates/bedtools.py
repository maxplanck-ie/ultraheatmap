import pandas as pd
import csv
import os
from pybedtools import BedTool

import gffutils

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

def extract_ge_folchange_per_peak(peaks, annotation, deseqtable, closestMapping):
    """

    """
    #### <START>
    ## The following set of functions produce a keymap (dictionary) mapping
    ## each gene to associated peaks (<1:n>-mapping, n > 0)
    ##

    ## keymap: peak_key:gene_id
    ## gene_id
    ## peak_keys: <chr>_<start>_<end>
    Peaks = BedTool(peaks)
    Peaks=Peaks.sort()    
    keyMap_closest = keymap_from_closest_genes(closestMapping, Peaks)
    peak2fc_table = extractFoldChange(keyMap_closest, deseqtable)

    return(peak2fc_table)


def extractFoldChange(keyMap_closest, deseqtable):
    """

    """
    newDF = pd.DataFrame()
    for key in keyMap_closest:
        seq = key.split('_')[0]
        start = key.split('_')[1]
        end = key.split('_')[2]
        geneExpression = deseqtable[deseqtable['GeneID'] == keyMap_closest[key]]
        geneExpression.insert(loc =0, column = 'peak_chr', value = seq)
        geneExpression.insert(loc =1, column = 'peak_start', value = start)
        geneExpression.insert(loc = 2, column = 'peak_end', value = end)
        
        newDF = newDF.append(geneExpression,ignore_index=True)
    return newDF


def parseGeneIdTable(table_file):
    return(pd.read_csv(table_file, sep ='\t'))
