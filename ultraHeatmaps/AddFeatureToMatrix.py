#!/usr/bin/env python

import os
import sys
import argparse
import yaml
import pandas as pd
import textwrap

sys.path.append(os.path.join(os.path.dirname(os.path.dirname((os.path.abspath(__file__))))))

from ultraHeatmaps.parseTables import extract_ge_folchange_per_peak, find_closest_genes, parseMatrixRegions, update_matrix_values, __read_tables_columns
from deeptools.heatmapper import heatmapper


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

def parse_args():
  """

  """

  desc = textwrap.dedent('''
    usage example:
    \tAddFeatureToMatrix -m <deepTools matrix> --output <combined deepTools matrix>
    \t\t--tables <tsv>
    \t\t--IDcolummn <name of id column>
    \t\t--Features <name of score column>

    The TSV is expected to contain header and minimum of two columns,
    an ID column and a value column.
    - ID column: name of region as in BED file used for compute[Ordered]Matrix
    - Score column: the score column that is added to the matrix

    If the closest gene has to be identified, an annotation GTF is expected,
    and the *ID column* needs to be the gene ID, as used in the annotation GTF.


    ''')
  parser=argparse.ArgumentParser(description = desc, formatter_class=CustomFormatter)

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
                      help="gene id tables or name based tables, tables should be space-separated.",
                      required = True)

  #optional arguments
  parser.add_argument("--Features",
                      "-f",
                      dest="Features",
                      nargs ='+',
                      help="A list of features of interest from gene id tables or name based tables",
                      default=["log2(FC)"])

  parser.add_argument("--IDcolumn",
                      dest="idcolumn",
                      type=str,
                      help="name of the column includes ids/names",
                      default="GeneID")


  parser.add_argument("--annotation",
                      "-a",
                      dest="Genome annotation file",
                      type=str,
                      metavar="STR",
                      default = None)
  parser.add_argument("--annotationFeature",
                      "-F",
                      dest="annotationFeature",
                      type=str,
                      help="annotation file can be filtered by a feature such as gene, exon or transcript",
                      default= None)
  parser.add_argument("--annotationOutput",
                      "-oa",
                      dest="annotationOutput",
                      type=str,
                      help="saving filtered annotation file if --annotationFeature",
                      default= None)

  parser.add_argument("--referencePoint",
                      dest="referencePoint",
                      type=str,
                      help="If closest TSS or TES is needed, otherwise closest gene body will be found",
                      default="TSS")

  return parser

def main():
   """
   Add a feature column to an existing deepTools matrix.

   The features are combined by the name column of the original BED file (use in
   computeMatrix).

   Alternatively, a genome annotation GTF can be provided, that allows to map each
   feature in the table to the closest gene.

   In each case, the original sorting of the deepTools matrix is kept.
   """

   parser = parse_args()
   args = parser.parse_args()

   #Check if the feature names are consistent between all the tables
   __read_tables_columns(args.tables,args.Features)

   hm = heatmapper()
   hm.read_matrix_file(args.deeptoolsMatrix)
   regions = parseMatrixRegions(hm.matrix.get_regions())
   #Using bedtool closest to map annotation and regions
   if args.annotation:
      closestMapping = find_closest_genes(regions, args.annotation, args.annotationFeature, args.annotationOutput,args.referencePoint) #XXX instead of all these arguments i can simply add args.

      # paste an extra column per table to the input matrix
      extract_ge_folchange_per_peak(regions, args.tables, closestMapping, args.Features, args.idcolumn,hm)


   else: #No closest gene is involved in this case , each enrichment id is individually checked and values are updated.
      update_matrix_values(regions, args.tables, args.Features, args.idcolumn,hm)
   #save the joint matrix obtained from either of cases
   hm.save_matrix(os.path.join(args.outputMatrix))

#if __name__ == "__main__":
#   main()
