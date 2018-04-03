import pandas as pd

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
    mapped=peaks.closest(sites, s=True, k=dictionary['distance']).saveas(output)


def __filter_annotation(dictionary): #XXX Here I would like to use a function from a package, if possible
   """
   filters annotation for the gene type of the interest
   """
   with open(dictionary['output']+"filtered.gtf","w") as filteredAnnotation:
        for feature in gffutils.DataIterator(dictionary['annotation']):
            if feature.featuretype == dictionary['featureToFilter']:
               filteredAnnotation.write(str(feature)+'\n')
          
   #anno = GffAnnotator(dictionary['annotation'], dictionary['featureToFilter'], "filtered", True)

def genes2Coordinates(geneids, gff_file, genome, version, fast_build = True):
    anno = GffAnnotator(gff_file, genome, version, fast_build)
    return(anno.geneId2Coordinates(geneids))

def parseGeneIdTable(table_file):
    return(pd.read_csv(table_file, sep ='\t'))

