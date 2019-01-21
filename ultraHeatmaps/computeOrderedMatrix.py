#!/usr/bin/env python

import os
import sys
import argparse
import yaml

sys.path.append(os.path.join(os.path.dirname(os.path.dirname((os.path.abspath(__file__))))))


#import necessary modules
import ultraHeatmaps.computeMatrix as cm

#Code directory:
configDir = os.path.dirname(os.path.realpath(__file__))
#Parse entered arguments
def parse_args(defaults={"kmeans":None, "hclust":None, "referencePoint":None, "metagene":None, "userconfig":None, "outFileSortedRegions": None, "refIndex" : None, "numberOfProcessors" : None}):

   """
   Parse the arguments from the command line
   """
   parser = argparse.ArgumentParser(description = "The program clusters regions and makes a matrix from the ordered regions.",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   # Required arguments
   parser.add_argument("-S",
                       "--Signal",
                       dest="bigwigs",
                       nargs='+',
                       help="bigwig files, the ordered matrix is computed"
                       "from.",
                       required=True)
   parser.add_argument("-R",
                       "--regions",
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

   #optional arguments
   parser.add_argument("--outFileSortedRegions",
                       dest="outFileSortedRegions",
                       help='[deepTools doc] File name in which the regions are '
                       'saved after skiping zeros or min/max threshold values. '
                       'The order of the regions in the file follows the sorting '
                       'order selected. This is useful, for example, to generate '
                       'other heatmaps keeping the sorting of the first heatmap. '
                       'Example: Heatmap1sortedRegions.bed',
                       default = defaults["outFileSortedRegions"])

   parser.add_argument("--kmeans",
                       dest="kmeans",
                       metavar="INT",
                       type=int,
                       help="number of clusters in k-means clustering",
                       default=defaults["kmeans"])
   parser.add_argument("--hclust",
                       dest="hclust",
                       metavar="INT",
                       type=int,
                       help="number of clusters in hierarchical clustering",
                       default=defaults["hclust"])

   parser.add_argument( '-i',
                        '--index',
                        "--referenceSampleIndex",
                       dest="refIndex",
                       nargs='+',
                       help='Index, 1-based, to define the reference samples. The '
                       'reference samples will be used for the clustering of the '
                       'matrix, before all samples are added. Space-separated',
                       default = defaults["refIndex"],
                       required=True)

   parser.add_argument("--metagene",
                       dest="metagene",
                       action="store_true",
                       help='[deepTools doc] When either a BED12 or GTF file are '
                       'used to provide regions, perform the computation on the '
                       'merged exons, rather than using the genomic interval '
                       'defined by the 5-prime and 3-prime most transcript bound '
                       '(i.e., columns 2 and 3 of a BED file). If a BED3 or BED6 '
                       'file is used as input, then columns 2 and 3 are used as an '
                        'exon.',
                       default=defaults["metagene"])
   parser.add_argument("-p",
                       "--numberOfProcessors",
                       dest="numberOfProcessors",
                       help='[deepTools doc] Number of processors to use. Type '
                       '"max/2" to use half the maximum number of processors or '
                       '"max" to use all available processors.',
                       type = int,
                       metavar="INT",
                       default=defaults["numberOfProcessors"])
   parser.add_argument("--referencePoint",
                       dest="referencePoint",
                       type=str,
                       help='[deepTools doc] The reference point for the plotting could be either'
                        'the region start (TSS), the region end (TES) or the'
                        'center of the region. Note that regardless of what you'
                        'specify, plotHeatmap/plotProfile will default to using'
                        '"TSS" as the label.',
                       default=defaults["referencePoint"])
   parser.add_argument("--config",
                       dest="userconfig",
                       help="Added to the default configuration, overwrites if "
                       "necessary.",
                       default=defaults["userconfig"])
   parser.add_argument('--samplesLabel',
                       help='[deepTools doc] Labels for the samples. This will then be passed to '
                       'plotHeatmap and plotProfile. The default is to use the '
                       'file name of the sample. The sample labels should be '
                       ' separated by spaces and quoted if a label itself'
                       'contains a space E.g. --samplesLabel label-1 "label 2"',
                        nargs='+')
   parser.add_argument('--cluster_mode',
                       dest = "cluster_mode",
                       help='The cluster is by default performed in the same way '
                       'as the final matrix same mode (reference-point/scale-regions). '
                       'To compute the clustering for \'scale-regions\', use \'0,0\' '
                       'for \'reference-point\', use \'A,B\', where A,B are up/downstream flanks [nt].',
                       type=str,
                       metavar="STR")

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
   if args.referencePoint:
       configfile["regionBodyLength"] = 0
       if configfile["beforeRegionStartLength"] == 0:
          configfile["beforeRegionStartLength"] = 1000
       if configfile["afterRegionStartLength"] == 0:
          configfile["afterRegionStartLength"] = 1000
   pre_cluster_mode =""
   boundries=[]
   if args.cluster_mode:
       a,b=args.cluster_mode.split(',')
       if b is '0' and a is '0':
           pre_cluster_mode = 'scale-regions'
       else:
           pre_cluster_mode = 'reference-point'
           boundries=[int(a),int(b)]
   add_diff(vars(args),configfile)
   #3. Generate an ordered region, using references only
   regions_list = args.regionOfInterest
   if args.refIndex:
       if configfile["outFileSortedRegions"] is None:
           path_name = os.path.dirname(os.path.abspath(args.matrixOutput))
           configfile["outFileSortedRegions"] = path_name+'/orderedBedFile.bed'
       cm.sortbyreference(args.regionOfInterest,args.refIndex,args.bigwigs,configfile, args, pre_cluster_mode, boundries)
       assert(os.path.getsize(configfile["outFileSortedRegions"]) > 0)
       regions_list = [configfile["outFileSortedRegions"]]

   #4.Build a matrix over all the samples
   hm = cm.computefinalmatrix(regions_list, args.bigwigs, configfile, args)

   matrix_output=os.path.join(args.matrixOutput)
   hm.save_matrix(matrix_output)


#if __name__ == "__main__":
    #main()
