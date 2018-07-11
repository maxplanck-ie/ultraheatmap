#!/usr/bin/env python

import os
import sys
import argparse
import yaml
import pandas as pd

sys.path.append(os.path.join(os.path.dirname(os.path.dirname((os.path.abspath(__file__))))))

from coordinates.parseTables import extract_ge_folchange_per_peak, find_closest_genes, parseMatrixRegions
from deeptools.heatmapper import heatmapper

def parse_args():
  """

  """
  parser=argparse.ArgumentParser()
  #required argumnets:
  parser.add_argument("--matrix",
                      "-m",
                      dest="deeptoolsMatrix",
                      type=str,
                      metavar="STR",
                      help="deeptools matrix",
                      required=True)
  parser.add_argument("--output",
                      "-om",
                      dest="outputMatrix",
                      type=str,
                      metavar="STR",
                      help="output matrix",
                      required=True)
  parser.add_argument("--tables",
                      "-t",
                      dest="tables",
                      nargs='+',
                      help="gene id tables or name based tables, names should be space-separated.",
                      required = True)

  #optional arguments
  parser.add_argument("--featureTypeToFilter",
                      "-F",
                      dest="featureToFilter",
                      type=str,
                      help="annotation file is filtered by gene, exon or transcriptt",
                      default= None)
  parser.add_argument("--annotationOutput",
                      "-oa",
                      dest="annotationOutput",
                      type=str,
                      help="filtered annotation file, it has to be added if featureTypeToFilter",
                      default= None)
  parser.add_argument("--annotation",
                      "-a",
                      dest="annotation",
                      type=str,
                      metavar="STR",
                      default = None)

  parser.add_argument("--Feature",
                      "-f",
                      dest="Feature",
                      type=str,
                      help="feature of interest from a gene id tables or name based tables",
                      default="log2(FC)")

  parser.add_argument("--IDcolumn",
                      dest="idcolumn",
                      type=str,
                      help="name of the column includes ids/names",
                      default="GeneID")

  return parser

def main():
   """
   Either the closest genes are foune and a deeptools-like matrix is created, if annotation file is provided,
   or a deeptools-like matrix directly from a provided enriched regions name-based files.
   In either case the output matrix is ordered and is appended to the input deeptools matrix.
   """
   parser = parse_args()
   args = parser.parse_args()
   hm = heatmapper()
   hm.read_matrix_file(args.deeptoolsMatrix)
   regions = parseMatrixRegions(hm.matrix.get_regions())
   #Using bedtool closest to map annotation and regions
   if args.annotation:
      closestMapping = find_closest_genes(regions, args.annotation, args.featureToFilter, args.annotationOutput)

      # paste an extra column per table to the input matrix
      extract_ge_folchange_per_peak(regions, args.tables, closestMapping, args.Feature, args.idcolumn,hm)


   #else:
       #TODO Thomas case, make the joined matrix

   #save the joint matrix obtained from either of cases
   hm.save_matrix(os.path.join(args.outputMatrix))

if __name__ == "__main__":
   main()
