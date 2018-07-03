#ComputeTFtrack.py

import os
import sys
import argparse
import yaml
import pandas as pd
import numpy as np
sys.path.append(os.path.join(os.path.dirname(os.path.dirname((os.path.abspath(__file__))))))

#import necessary modules
from pybedtools import BedTool

from deeptoolsapi.computeMatrix import compute_matrix
from deeptoolsapi.deeptoolsMatrix import Matrix, read_matrix_file
from coordinates.bedtools import find_closest_genes, __parseRegions, __getValuesFromDEseqTable
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

  parser.add_argument("--transcriptionFactors", #comma separated TF files
                      "-tf",
                      dest="transcriptionFactors",
                      type=str,
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
                      default="percentage")


  return parser

def main():
    """
    For each TF file one column is generated to be added to the deeptools matrix
    """
    parser = parse_args()
    args = parser.parse_args()
    TFfiles = args.transcriptionFactors
    files_list =[str(filename) for filename in TFfiles.split(',')]
    count = 0
    regionFile = pd.read_csv(args.regionOfInterest, sep ='\t')
    regions=[]
    names = []
    for id, row in regionFile.iterrows():
         string= str(row[0])+";"+str(row[1])+";"+str(row[2])+";"+str(row[0])+":"+str(row[1])+"-"+str(row[2])+";.;.;"
         regions.append(string)
         names.append(str(row[0])+":"+str(row[1])+"-"+str(row[2]))
    
    valuesTab = np.empty((len(regions), len(files_list)), dtype=float)
    for i, table in enumerate(files_list):
        tf_score = pd.read_csv(table,sep = '\t')
        values = __getValuesFromDEseqTable(names, tf_score, args.Feature)
        valuesTab[:,i] = values
    matrix_output = os.path.join(args.output, "TF.matrix.gz")
    matrix = Matrix(regions = __parseRegions(regions) , matrix = valuesTab, group_boundaries = [0,len(regions)], sample_boundaries =  [x for x in range(0, len(files_list) + 1, 1)])
    matrix.save_matrix(matrix_output)

if __name__ == "__main__":
   main()
