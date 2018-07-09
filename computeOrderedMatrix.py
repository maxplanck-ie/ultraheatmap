#!/usr/bin/env python

import os
import sys
import argparse
import yaml

sys.path.append(os.path.join(os.path.dirname(os.path.dirname((os.path.abspath(__file__))))))


#import necessary modules
import deeptoolsapi.computeMatrix as cm

#Code directory:
code_dir = os.path.dirname(os.path.realpath(__file__))

#Parse entered arguments
def parse_args(defaults={"kmeans":None, "hclust":None, "referencePoint":None, "metagene":None, "userconfig":None, "outFileSortedRegions": None, "refIndex" : None}):

   """
   Parse the arguments from the command line
   """
   parser = argparse.ArgumentParser()
   # Required arguments
   parser.add_argument("-b",
                       "--Signal",
                       dest="bigwigs",
                       nargs='+',
                       help="All the Bigwig files in spaced format",
                       required=True)
   parser.add_argument("-R",
                       "--regions",
                       dest="regionOfInterest",
                       nargs='+',
                       help="All regions of interest, spaced .bed/.gtf files",
                       required=True)
   parser.add_argument("-om",
                    "--matrixOutput",
                    dest="matrixOutput",
                    help="the output matrix",
                    required=True)

   #optional arguments
   parser.add_argument("-os",
                       "--sortedOutput",
                       dest="outFileSortedRegions",
                       help="the ordered region file (.bed)",
                       default = defaults["outFileSortedRegions"])
   parser.add_argument("--kmeans",
                       dest="kmeans",
                       metavar="INT",
                       type=int,
                       help="number of k-means clusters",
                       default=defaults["kmeans"])
   parser.add_argument("--hclust",
                       dest="hclust",
                       metavar="INT",
                       type=int,
                       help="number of clusters in hierarchical clustering",
                       default=defaults["hclust"])
   parser.add_argument("--refIndex",
                       dest="refIndex",
                       nargs='+',
                       help="Indices of bigwig files which pointed to the references. Several indices can be separated by space. Note that the numbers are zero based!",
                       default = defaults["refIndex"])
   parser.add_argument("--metagene",
                       dest="metagene",
                       action="store_true",
                       help="when region is .GTF or .BED12 and mode is scale-regions",
                       default=defaults["metagene"])
   parser.add_argument("--refPoint",
                       dest="referencePoint",
                       type=str,
                       help="Reference point for plotting can be set to TSS, TES or center",
                       default=defaults["referencePoint"])
   parser.add_argument("--config",
                       dest="userconfig",
                       help="this will be added to the dieafult config file",
                       default=defaults["userconfig"])

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
   with open(os.path.join(code_dir, "computeOrderedMatrix.yaml"), 'r') as stream:
     defaultconfigfile = yaml.load(stream)
   #2. Parse the arguments
   parser = parse_args(defaultconfigfile)
   args = parser.parse_args()
   #2.1 modify config file if needed
   configfile=defaultconfigfile
   add_diff(vars(args),configfile)
   if args.userconfig:
      configfile= merge_dictionaries(configfile, args.userconfig)

   #3. Generate an ordered region, using references only
   if args.refIndex:
       cm.sortbyreference(args.regionOfInterest,args.refIndex,args.bigwigs,configfile)
       if os.path.getsize(configfile["outFileSortedRegions"]) > 0:
          regions_list = [configfile["outFileSortedRegions"]]

   #4.Build a matrix over all the samples
   hm = cm.computefinalmatrix(regions_list, args.bigwigs, configfile)

   matrix_output=os.path.join(args.matrixOutput)
   hm.save_matrix(matrix_output)








if __name__ == "__main__":
    main()
