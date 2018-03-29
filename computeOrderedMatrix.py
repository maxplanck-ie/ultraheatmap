#!/usr/bin/env python

import os
import sys
import argparse
import yaml

sys.path.append(os.path.join(os.path.dirname(os.path.dirname((os.path.realpath(__file__))))))


#import necessary modules
import deeptoolsapi.computeMatrix as cm

#Code directory:
code_dir = os.path.dirname(os.path.realpath(__file__))

#Parse entered arguments
def parse_args(defaults={"mode":None, "flanking_str":None, "unscaled_str":None, "kmeans_clust":None, "hclust":None, "refPoint":None, "extramatrix":None, "metagene":None, "userconfig":None}):

   """
   Parse the arguments from the command line
   """
   parser = argparse.ArgumentParser()
   # Required arguments
   parser.add_argument("-b",
                       "--Signal",
                       dest="bigwigs",
                       type=str,
                       help="All the Bigwig files in comma separated format",
                       required=True)
   parser.add_argument("-R",
                       "--regions",
                       dest="regionOfInterest",
                       help="Region of interest(.bed/.gtf)",
                       required=True)
   parser.add_argument("-o",
                       "--output",
                       dest="outputDir",
                       help="the directory that all the output will besaved",
                       required=True)
   #optional arguments
   parser.add_argument("-m",
                       "--mode",
                       dest="mode",
                       type = str,
                       help= "available modes are reference-point and scale-regions",
                       default=defaults["mode"])
   parser.add_argument("-f",
                      "--flank",
                      dest="flanking_str",
                      metavar="INT",
                      type=int,
                      help="flanking regions",
                      default=defaults["flanking_str"])
   parser.add_argument("--unscaled_str",
                       dest="unscaled_str",
                       type=int,
                       metavar="INT",
                       help="number of bases from 5 prime and 3 prime to exclude from the region",
                       default=defaults["unscaled_str"])
   parser.add_argument("--kmeans_clust",
                       dest="kmeans_clust",
                       metavar="INT",
                       type=int,
                       help="number of k-means clusters",
                       default=defaults["kmeans_clust"])
   parser.add_argument("--hclust",
                       dest="hclust",
                       metavar="INT",
                       type=int,
                       help="number of clusters in hierarchical clustering",
                       default=defaults["hclust"])
   parser.add_argument("--refIndex", #XXX Note that it is zero based!
                       dest="refIndex",
                       type=str,
                       metavar="STRING",
                       help="bigwig files pointed to references in comma separated format")
   parser.add_argument("--metagene",
                       dest="metagene",
                       action="store_true",
                       help="when region is .GTF or .BED12 and mode is scale-regions",
                       default=defaults["metagene"])
   parser.add_argument("--refPoint",
                       dest="refPoint",
                       type=str,
                       help="Reference point for plotting can be set to TSS, TES or center",
                       default=defaults["refPoint"])
   parser.add_argument("--addmatrix",
                       "-M",
                       dest="extramatrix",
                       type=str,
                       metavar="STRING",
                       help="add a matrix",
                       default=defaults["extramatrix"])
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
   Main function to integrate heatmaps (and to add external tracks??)
   """
   # First part of the code is applied for all cases
   defaultconfigfile = {}
   #1. Read the config file
   with open(os.path.join(code_dir, "computeOrderedMatrix.yaml"), 'r') as stream:
     defaultconfigfile = yaml.load(stream)
   #2. Parse the arguments
   parser = parse_args(defaultconfigfile)
   args = parser.parse_args()
   output_dir =os.path.abspath(args.outputDir)
   # modify config file if needed
   configfile=defaultconfigfile
   add_diff(vars(args),configfile)
   print(configfile)
   if args.userconfig:
      configfile= merge_dictionaries(configfile, args.userconfig)
   with open(os.path.join(output_dir,'configfile.yaml'), 'w') as c:
      yaml.dump(configfile, c, default_flow_style=False)

   regionOfInterest =os.path.abspath(args.regionOfInterest)
   open(regionOfInterest,'r')

   #3. Generate an ordered region, using references onyl
   bigwig_list=[str(bw) for bw in args.bigwigs.split(',')] ##This can later on be used for refIndex
   bigwig_files = " ".join(bigwig_list)

   regions = regionOfInterest
   if args.refIndex:
       orderedbed = cm.sortbyreference(regions,args.refIndex,bigwig_list,configfile)
       regions = orderedbed
       assert regions == os.path.join(output_dir,"ordered.bed")

   #4.Built matrices over all the samples, add closest gene matrix is provided

   cm.computefinalmatrix(regions, bigwig_files, configfile)









if __name__ == "__main__":
    main()
