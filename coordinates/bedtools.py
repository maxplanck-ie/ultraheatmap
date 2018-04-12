import pandas as pd

import pybedtools

from gffannotator.gffannotator import GffAnnotator

def find_closest_genes(dictionary,output):
    """
    Find the closest gene using bedtools.closest
    """
    peaks=pybedtools.BedTool(dictionary['regionOfInterest']).sort().saveas()
    sites=pybedtools.BedTool(dictionary['annotation']).sort().saveas()
    if dictionary['geneToFilter'] != None:
       print(dictionary['geneToFilter'])
       __filter_annotation(dictionary)
       sites=pybedtools.BedTool(dictionary['output']+"filtered.gtf").sort().saveas()
    mapped=sites.closest(peaks, s=True, k=1000).saveas(output) #TODO k should become an argument


def __filter_annotation(dictionary): #XXX Here I would like to use a function from a package, if possible
   """
   filters annotation for the gene type of the interest
   """
   with open(dictionary['output']+"filtered.gtf","w") as filteredAnnotation:
        with open(dictionary['annotation'], "r") as f:
             line = f.readline()
             while line:
                if dictionary['geneToFilter'] in line: #TODO this is wrong! need to read recorde by record rather than line by line!
                   filteredAnnotation.write(line)
                line = f.readline()

def genes2Coordinates(geneids, gff_file, genome, version, fast_build = True):
    anno = GffAnnotator(gff_file, genome, version, fast_build)
    return(anno.geneId2Coordinates(geneids))

def parseGeneIdTable(table_file):
    return(pd.read_csv(table_file, sep ='\t'))
