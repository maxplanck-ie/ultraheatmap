#!/usr/bin/env python

import os
import sys
import argparse
import yaml
import tempfile

sys.path.append(os.path.join(os.path.dirname(os.path.dirname((os.path.abspath(__file__))))))


#import necessary modules
import ultraheatmap.computeMatrix as cm

#Code directory:
configDir = os.path.dirname(os.path.realpath(__file__))
#Parse entered arguments
def parse_args(defaults={}):

   """
   Parse the arguments from the command line
   """
   parser = argparse.ArgumentParser(description = "The program clusters regions and makes a matrix from the ordered regions.",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   # Required arguments
   parser.add_argument("-S",
                       "--scoreFileName",
                       dest="bigwigs",
                       nargs='+',
                       help="bigwig files, the ordered matrix is computed"
                       "from.",required=True)
   parser.add_argument("-R",
                       "--regionsFileName",
                       dest="regionOfInterest",
                       nargs='+',
                       help="BED files definig the genomic regions of the matrix."
                       "Multiple files can be provided, but the per group information"
                       "will be lost due to the clustering",
                       required=True)
   parser.add_argument("-o",
                       "--outFileName",
                       dest="matrixOutput",
                       help="Matrix clustered by the given reference samples",
                       required=True)

   parser.add_argument( "-i",
                        "--index",
                        "--referenceSampleIndex",
                        dest="refIndex",
                        nargs='+',
                        help='Index, 1-based, to define the reference samples. The '
                        'reference samples will be used for the clustering of the '
                        'matrix, before all samples are added. Space-separated',
                        required=True)

   #optional arguments
   parser.add_argument("-p",
                       "--numberOfProcessors",
                       dest="numberOfProcessors",
                       help='[deepTools doc] Number of processors to use. Type '
                       '"max/2" to use half the maximum number of processors or '
                       '"max" to use all available processors.',
                       type = int,
                       metavar="INT",
                       default=1)

   parser.add_argument("--outFileSortedRegions",
                       dest="outFileSortedRegions",
                       help='[deepTools doc] File name in which the regions are '
                       'saved after skiping zeros or min/max threshold values. '
                       'The order of the regions in the file follows the sorting '
                       'order selected. This is useful, for example, to generate '
                       'other heatmaps keeping the sorting of the first heatmap. '
                       'Example: Heatmap1sortedRegions.bed',
                       default = None)
   parser.add_argument("--outputReferenceMatrix",
                        dest="outputReferenceMatrix",
                        help="Matrix on the reference sampels only before clustering",
                        default=None)

   parser.add_argument("--kmeans",
                       dest="kmeans",
                       metavar="INT",
                       type=int,
                       help="number of clusters in k-means clustering",
                       default=None)

   parser.add_argument("--hclust",
                       dest="hclust",
                       metavar="INT",
                       type=int,
                       help="Number of clusters to compute using hierarchical clustering as defined by deepTools plotHeatmap",
default=None)


   parser.add_argument("--config",
                       dest="userconfig",
                       help="Added to the default configuration, overwrites if "
                       "necessary.",
                       default=None)

   return  parser



def merge_dictionaries(a,b):
    """
    merging two dictionaries
    """
    merged = {**a, **b}
    return merged


def add_diff(a,b):
   """
   Add the difference between the default config file and argumnets to the configfile. Order matters!
   b is defaultconfig
   a is vars(args)
   """
   for key in b.keys():
       if key in a:
          b[key]=a[key]


def main():
   """
   Compute a matrix from an ordered region file
   """
   # First part of the code is applied for all cases
   defaultconfigfile = {}
   #1. Read the config file
   with open(os.path.join(configDir, 'configs' ,"computeOrderedMatrix.yaml"), 'r') as stream:
     defaultconfigfile = yaml.load(stream)
   #2. Parse the arguments
   parser = parse_args(defaultconfigfile)
   args = parser.parse_args()
   #2.1 modify config file if needed
   configfile=defaultconfigfile
   if args.userconfig:
       with open(os.path.join(args.userconfig), 'r') as stream:
            userconfigfile = yaml.load(stream)
            configfile= merge_dictionaries(configfile, userconfigfile)

   configfile= merge_dictionaries(configfile, vars(args))
   configfile['numberOfProcessors'] = args.numberOfProcessors

   #3. Generate an ordered region, using references only
   if configfile["outFileSortedRegions"] is None:
       path_name = os.path.dirname(os.path.abspath(args.matrixOutput))
       configfile["outFileSortedRegions"] = path_name+'/orderedBedFile.bed'


   cm.sortbyreference(configfile["regionOfInterest"], configfile["bigwigs"], configfile["refIndex"], configfile)
   assert(os.path.getsize(configfile["outFileSortedRegions"]) > 0)

   #4.Build a matrix over all the samples
   hm = cm.computefinalmatrix(configfile["outFileSortedRegions"], configfile["bigwigs"], configfile)

   matrix_output=os.path.join(args.matrixOutput)
   hm.save_matrix(matrix_output)
