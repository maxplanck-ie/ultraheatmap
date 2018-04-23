import pandas as pd
import csv
from pybedtools import BedTool

import gffutils

from coordinates.mapClosestGenes import keymap_from_closest_genes

def find_closest_genes(peaks, annotation, dictionary, filename = None):
    """
    Find the closest gene using bedtools.closest
    """
    peaks=peaks.sort()
    sites=annotation.sort()
    if dictionary['featureToFilter'] != None:
        print(dictionary['featureToFilter'])
        __filter_annotation(dictionary)
        sites=BedTool(dictionary['output']+"filtered.gtf").sort()

    mapped=peaks.closest(sites, t="first")

    if filename:
        mapped.saveas(filename)

    return(mapped)

def __filter_annotation(dictionary):
   """
   filters annotation for the gene type of the interest
   """
   with open(dictionary['output']+"filtered.gtf","w") as filteredAnnotation:
        for feature in gffutils.DataIterator(dictionary['annotation']):
            if feature.featuretype == dictionary['featureToFilter']:
               filteredAnnotation.write(str(feature)+'\n')

def extract_ge_folchange_per_peak(peaks, annotation, deseqtable):
    """

    """
    #### <START>
    ## The following set of functions produce a keymap (dictionary) mapping
    ## each gene to associated peaks (<1:n>-mapping, n > 0)
    ##
    ## keymap: peak_key:gene_id
    ## gene_id
    ## peak_keys: <chr>_<start>_<end>

    peaks = BedTool(dictionary['regionsOfInterest'])
    annotation = BedTool(dictionary['annotation'])
    closestMapping = find_closest_genes(peaks, annotation, dictionary, None)
    keyMap_closest = keymap_from_closest_genes(closestMapping, peaks)

    # TODO
    peak2fc_table = extractFoldChange(keyMap_closest, deseqtable)

    return(peak2fc_table)


def parseGeneIdTable(table_file):
    return(pd.read_csv(table_file, sep ='\t'))
