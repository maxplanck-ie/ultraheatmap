import pandas as pd
import csv
import os
import sys
from pybedtools import BedTool
import numpy as np
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

def extract_ge_folchange_per_peak(peaks, annotation, deseqtables, closestMapping,deseqfeature):
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
    peak2fc_table = extractFoldChange(keyMap_closest, deseqtables,deseqfeature)

    return(peak2fc_table)


def extractFoldChange(keyMap_closest, deseqtables,deseqfeature):
    """

    """
    geneIdtables=[]
    for table in deseqtables:
       geneIdtables.append(pd.read_csv(table,sep ='\t', squeeze = True))
    newDF = pd.DataFrame()
    for key in keyMap_closest:
        seq,start,end = key.split('_')[:3]
        geneExpression=pd.DataFrame()
        for i, table in enumerate(geneIdtables):
            colname=deseqfeature+"_deseq"+str(i)
            if keyMap_closest[key] in table['GeneID'].values:
               geneExpression.insert(loc = i, column = colname, value = table[table['GeneID'] == keyMap_closest[key]][deseqfeature])
            else:
               geneExpression.insert(loc = i, column = colname, value = np.nan)
        geneExpression.insert(loc =0, column = 'peak_chr', value = seq)
        geneExpression.insert(loc =1, column = 'peak_start', value = start)
        geneExpression.insert(loc = 2, column = 'peak_end', value = end)
        newDF = newDF.append(geneExpression,ignore_index=True)
    return newDF


def parseGeneIdTable(table_file):
    return(pd.read_csv(table_file, sep ='\t'))
