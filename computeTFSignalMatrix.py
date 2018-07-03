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
    regions = pd.read_csv(args.regionOfInterest, sep ='\t')
    values_table = np.empty((len(regions), len(files_list)), dtype=float)
    for file in files_list:
        tf_score = pd.read_csv(file,sep = '\t')
        values =[]
        regs = []
        for id, row in regions.iterrows():
             string= str(row[0])+":"+str(row[1])+"-"+str(row[2])
             if string in tf_score['name'].values:
                 x = float(tf_score[tf_score.name == string][args.Feature])
                 if np.isnan(x):
                    x = np.nan
                 values += [ x ]
             else:
                 values += [ np.nan ]
             strand = "."
             score = "."
             reg = [(int(row[1]),int(row[2]))]
             regs.append([row[0], reg, string, len(regions), strand, score])
        values_table[:,count] = values
        count +=1
    matrix_output = os.path.join(args.output, "TF.matrix.gz")
    matrix = Matrix(regs , values_table, group_boundaries = [0,len(regions)], sample_boundaries =  [x for x in range(0, len(files_list) + 1, 1)])
    matrix.save_matrix(matrix_output)

if __name__ == "__main__":
   main()
