#ComputeTFtrack.py

import os
import sys
import argparse
import yaml
import pandas as pd
import numpy as np
from itertools import accumulate
sys.path.append(os.path.join(os.path.dirname(os.path.dirname((os.path.abspath(__file__))))))

#import necessary modules
from pybedtools import BedTool

from deeptoolsapi.computeMatrix import compute_matrix
from coordinates.bedtools import find_closest_genes, __parseRegions, __getValuesFromDEseqTable
from coordinates.bedtools import extract_ge_folchange_per_peak

def parse_args():
  """

  """
  parser=argparse.ArgumentParser()
  #required argumnets:
  parser.add_argument("--region", #comma separated region files - regions cannot have header!! XXX Always a uniques name is needed
                      "-R",
                      dest="regionOfInterest",
                      type=str,
                      metavar="STR",
                      help="Regions of intereset(.bed/.gtf)",
                      required=True)

  parser.add_argument("--transcriptionFactors", #comma separated TF files
                      "-tf",
                      dest="transcriptionFactors",
                      type=str,
                      metavar="STR",
                      help="transcription factor signals (.bed)",
                      required = True)

  parser.add_argument("--output",
                      "-o",
                      dest="output",
                      type=str,
                      metavar="STR",
                      required=True)


  parser.add_argument("--Feature",
                      dest="Feature",
                      type=str,
                      help="feature of interest",
                      default="ratio")

  parser.add_argument("--IDcolumn",
                      dest="idcolumn",
                      type=str,
                      help="name of the column includes ids/names",
                      default="GeneID")

  return parser

def main():
    """
    For each TF file one column is generated to be added to the deeptools matrix
    """
    parser = parse_args()
    args = parser.parse_args()
    TFfiles = args.transcriptionFactors
    files_list =[str(filename) for filename in TFfiles.split(',')]
    groupBoundries = [0]
    names = []
    regions=[]
    regions_list =[str(filename) for filename in args.regionOfInterest.split(',')]
    for region in regions_list:
        regionFile = pd.read_csv(region, sep ='\t', header=None) ##XXX IS it always without header?
        rows, columns = regionFile.shape
        assert columns >= 3
        assert columns < 7
        if columns == 3:
            regionFile[len(regionFile.columns)] = "."
            rows, columns = regionFile.shape
        if columns == 4:
           regionFile[len(regionFile.columns)] = "."
           rows, columns = regionFile.shape
        if columns == 5:
            regionFile[len(regionFile.columns)] = "."
            rows, columns = regionFile.shape
        for id, row in regionFile.iterrows():
            if row[3]  == ".":
               row[3] = str(row[0])+":"+str(row[1])+"-"+str(row[2])
            string= str(row[0])+";"+str(row[1])+";"+str(row[2])+";"+str(row[3])+";"+str(row[4])+";"+str(row[5])+";"
            regions.append(string)
            names.append(row[3])
        groupBoundries.append(len(regionFile))

    valuesTab = np.empty((len(regions), len(files_list)), dtype=float)
    print(valuesTab.shape)
    for i, table in enumerate(files_list):
       tf_score = pd.read_csv(table,sep = '\t')
       values = __getValuesFromTable(names, tf_score, args.Feature, args.idcolumn)
       print(len(values))
       valuesTab[:,i] = values

    matrix_output = os.path.join(args.output, "TF.matrix.gz")
    matrix = Matrix(regions = __parseRegions(regions) , matrix = valuesTab, group_boundaries = list(accumulate(groupBoundries)), sample_boundaries =  [x for x in range(0, len(files_list) + 1, 1)])
    matrix.save_matrix(matrix_output)

if __name__ == "__main__":
   main()
