import pandas as pd
import csv
import os
from pybedtools import BedTool

import gffutils

from gffannotator.gffannotator import GffAnnotator
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

def map_peaks_to_geneID(dictionary):
    """

    """
    #### <START>
    ## The following set of functions produce a keymap (dictionary) mapping
    ## each gene to associated peaks (<1:n>-mapping, n > 0)
    ##
    ## keymap keys <gene key>:[<peak keys>] (multimap)
    ## gene keys: gff gene_id
    ## peak keys: <chr>_<start>_<end>

    peaks = BedTool(dictionary['regionOfInterest'])
    annotation = BedTool(dictionary['annotation'])
    closestMapping = find_closest_genes(peaks, annotation, dictionary)
    
    keyMap_closest = keymap_from_closest_genes(closestMapping, peaks)
    #XXX It seems it would save us so much time if the featuredb is built over a gff file that contains only the closest genes. Doest it make sense?
    geneID2coordMap = genes2Coordinates(dictionary['geneIDtable'], dictionary['annotation'])
    print(geneID2coordMap)

def genes2Coordinates(geneids, gff_file, fast_build = True):
    anno = GffAnnotator(gff_file, fast_build)
    
    ids = open(geneids, 'r')
    next(ids)
    lines = ids.readlines()
    id_list=list()
    for l in lines:       
        id_list.append(l.split()[0])  
    return (anno.geneid2keymap(id_list))


def parseGeneIdTable(table_file):
    return(pd.read_csv(table_file, sep ='\t'))
