#!/usr/bin/env python

import os
import sys
import argparse
import yaml

sys.path.append(os.path.join(os.path.dirname(os.path.dirname((os.path.realpath(__file__))))))

#import necessary modules
from deeptoolsapi.computeMatrix import compute_matrix
from coordinates.bedtools import find_closest_genes
from coordinates.bedtools import map_peaks_to_geneID

def parse_args():
  """

  """
  parser=argparse.ArgumentParser()
  #required argumnets:
  parser.add_argument("--region",
                      "-R",
                      dest="regionOfInterest",
                      type=str,
                      metavar="STR",
                      help="Regions of intereset(.bed/.gtf)",
                      required=True)
  parser.add_argument("--signal",
                      "-b",
                      dest="bigwigs",
                      type=str,
                      metavar="STR",
                      help="All the Bigwig files in comma separated format",
                      required=True)
  parser.add_argument("--annotation",
                      "-a",
                      dest="annotation",
                      type=str,
                      metavar="STR",
                      required=True)
  parser.add_argument("--output",
                      "-o",
                      dest="output",
                      type=str,
                      metavar="STR",
                      required=True)
  parser.add_argument("--mode",
                      "-m",
                      dest="mode",
                      type=str,
                      metavar="STR",
                      help="scale-regions or reference-point")
  parser.add_argument("-f",
                      "--flank",
                      dest="flanking_str",
                      metavar="INT",
                      type=int,
                      help="flanking regions",
                      default=1000)
  parser.add_argument("--unscaled_str",
                      dest="unscaled_str",
                      type=int,
                      metavar="INT",
                      help="number of bases from 5 prime and 3 prime to exclude from the region",
                      default=100)
  parser.add_argument("--kmeans_clust",
                      dest="kmeans_clust",
                      metavar="INT",
                      type=int,
                      help="number of k-means clusters",
                      default=0)
  parser.add_argument("--hclust",
                      dest="hclust",
                      metavar="INT",
                      type=int,
                      help="number of clusters in hierarchical clustering",
                      default=0)
  parser.add_argument("--metagene",
                      dest="metagene",
                      action="store_true",
                      help="when region is .GTF or .BED12 and mode is scale-regions",
                      default=False)
  parser.add_argument("--refPoint",
                      dest="refPoint",
                      type=str,
                      help="Reference point for plotting can be set to TSS, TES or center",
                      default='center')
  parser.add_argument("--featureTypeToFilter",
                      "-F",
                      dest="featureToFilter",
                      type=str,
                      help="annotation file is filtered by gene, exon or transcriptt",
                      default=None)
  parser.add_argument("--geneIDs",
                      "-g",
                      dest="geneIDtable",
                      type=str,
                      help="gene id table",
                      default=None)

  return parser

def main():
   """
   Main function to find the closest genes and to generate a matrix
   """
   parser = parse_args()
   args = parser.parse_args()
   print(args)
   #Using bedtool closest to map annotation and regions
   map_peaks_to_geneID(vars(args))
   #find_closest_genes(vars(args),args.output+"closestGene.bed")
   #compute_matrix is run over mapped.bed and .bw files
   if args.mode:
      bigwig_list=[str(bw) for bw in args.bigwigs.split(',')]
      bigwig_files = " ".join(bigwig_list)
      matrix_output=os.path.join(args.output, args.mode+"_nearest_gene.matrix")

      compute_matrix( vars(args)['mode'], bigwig_files, args.output+"mapped.bed",matrix_output, vars(args))


if __name__ == "__main__":
    main()
