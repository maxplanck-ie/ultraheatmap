#!/usr/bin/env python

import os
import sys
import argparse
import yaml
import pandas as pd

sys.path.append(os.path.join(os.path.dirname(os.path.dirname((os.path.realpath(__file__))))))

#import necessary modules
from pybedtools import BedTool

from deeptoolsapi.computeMatrix import compute_matrix
from deeptoolsapi.deeptoolsMatrix import Matrix, read_matrix_file
from coordinates.bedtools import find_closest_genes
from coordinates.bedtools import extract_ge_folchange_per_peak

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

  parser.add_argument("--geneIDs", #comma separated gene id tables
                      "-g",
                      dest="geneIDtables",
                      type=str,
                      help="gene id table",
                      required = True)

  #optional arguments
  parser.add_argument("--featureTypeToFilter",
                      "-F",
                      dest="featureToFilter",
                      type=str,
                      help="annotation file is filtered by gene, exon or transcriptt",
                      default='exon')

  parser.add_argument("--deseqFeature",
                      dest="deseqFeature",
                      type=str,
                      help="feature of interest from a deseq table",
                      default="log2(FC)")


  return parser

def main():
   """
   Main function to find the closest genes and to generate a matrix
   """
   parser = parse_args()
   args = parser.parse_args()
   geneIDtables = args.geneIDtables
   table_list =[str(geneId) for geneId in geneIDtables.split(',')]

   #Using bedtool closest to map annotation and regions
   closestMapping = find_closest_genes(args.regionOfInterest, args.annotation, args.featureToFilter,args.output)

   # deeptoolsMatrix(peak2foldchange)
   peak2foldchange = extract_ge_folchange_per_peak(args.regionOfInterest, args.annotation, table_list, closestMapping, args.deseqFeature)
   matrix_output=os.path.join(args.output, "closestGene.matrix.gz")
   peak2foldchange.save_matrix(matrix_output)

if __name__ == "__main__":
   main()
