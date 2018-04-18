import pandas as pd
import csv
import pybedtools

import gffutils

from gffannotator.gffannotator import GffAnnotator

def find_closest_genes(dictionary,output):
    """
    Find the closest gene using bedtools.closest
    """
    peaks=pybedtools.BedTool(dictionary['regionOfInterest']).sort().saveas()
    sites=pybedtools.BedTool(dictionary['annotation']).sort().saveas()
    if dictionary['featureToFilter'] != None:
       print(dictionary['featureToFilter'])
       __filter_annotation(dictionary)
       sites=pybedtools.BedTool(dictionary['output']+"filtered.gtf").sort().saveas()
    mapped=peaks.closest(sites, t="first").saveas(output)


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

    genes2Coordinates(dictionary['geneIDtable'], dictionary['annotation'])


def genes2Coordinates(geneids, gff_file, fast_build = True):
    anno = GffAnnotator(gff_file, fast_build)
    ids = open(geneids, 'r')
    next(ids)
    lines = ids.readlines()
    id_list=list()
    for l in lines:
        id_list.append(l.split()[0])
    mappedAnno=anno.geneId2Coordinates(id_list) #TODO
    


def parseGeneIdTable(table_file):
    return(pd.read_csv(table_file, sep ='\t'))
