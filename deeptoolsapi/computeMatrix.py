import subprocess
import shlex
import sys
import os
import pandas as pd
import numpy as np

sys.path.append(os.path.join(os.path.dirname(os.path.dirname((os.path.realpath(__file__))))))
from bashwrapper.bashwrapper import Bash
from deeptoolsapi.plotHeatMap import plot_heatmap
from deeptoolsapi.deeptoolsMatrix import read_matrix_file
#from deeptools.heatmapper import heatmapper as dh #TODO
#Loading deeptools:
deeptools_module_load="module load deeptools"
def __computeScaledRegions(bigwigs, regions, matrix, configfile):
    opt_metagene=""
    if configfile['metagene']:
        opt_metagene = '--metagene'
    cmd = "computeMatrix scale-regions -S {bigwigs} -R {regions} -a {flank} -b {flank} --unscaled5prime {unscaled5} --unscaled3prime {unscaled3} {metagene} -o {outputmatrix}".format(
                bigwigs= bigwigs,
                regions=' '.join(regions.split(',')),
                flank = str(configfile['flanking_str']),
                unscaled5 = str(configfile['unscaled_str']),
                unscaled3 = str(configfile['unscaled_str']),
                metagene = opt_metagene,
                outputmatrix = matrix)
    return(cmd)
#     dh.computeMatrix()
def __computeReferencePoint(bigwigs, regions, matrix, configfile):
    opt_metagene=""
    if configfile['metagene']:
        opt_metagene = '--metagene'
    cmd = "computeMatrix reference-point --referencePoint {refPoint} -S {bigwigs} -R {regions} -a {flank} -b {flank} {metagene} -o {outputmatrix}".format(
                refPoint=configfile['refPoint'],
                bigwigs=' '.join(bigwigs.split(',')),
                regions=' '.join(regions.split(',')),
                flank = str(configfile['flanking_str']),
                metagene = opt_metagene,
                outputmatrix = matrix)
    return(cmd)

#first a matrix is computed using deeptools/computeMatrix, then deeptools/cbind is used to merge the generated matrix with the closest-genes matrix if the closest-genes matrix is provided.
def compute_matrix(mode, bw, bed, matrix, configfile):
   """
   computing the corresponding matrix using deeptools/computeMatrix
   """
   computeMatrix_cmd = [deeptools_module_load]
   if mode == "scale-regions":
       computeMatrix_cmd.append( __computeScaledRegions(bw, bed, matrix, configfile) )
#       __computeScaledRegions(bw, bed, matrix, configfile)
   else:
       assert mode == "reference-point"
       computeMatrix_cmd.append( __computeReferencePoint(bw, bed, matrix, configfile) )


   cmd=";".join(computeMatrix_cmd)
#  print("cmd\n"+str(cmd))
   computeMatrix_bash = Bash(cmd)


def sortbyreference(regions,refIndex,bigwig_list,configfile):
    indexList=[int(index) for index in refIndex.split(',')]
    refList=[]
    for index in indexList:
        refList.append(bigwig_list[index])
    #only on the ref.ones
    ref_bw=" ".join(refList)
    matrix_name = 'refonly_matrix'
    matrix_output=os.path.join(configfile['outputDir'], matrix_name)
    compute_matrix(configfile['mode'], ref_bw, regions,matrix_output, configfile)
    plot_heatmap(configfile['mode'], matrix_output, configfile,True)
    orderedbed=os.path.join(configfile['outputDir'],"ordered.bed")
    os.path.exists(os.path.join(configfile['outputDir'],"ordered.bed"))
    return orderedbed


def __cbind_matrix(matrix1, matrix2, output_dir):
    cbind_cmd = [deeptools_module_load]
    cbind_cmd.append("computeMatrixOperations cbind -m {matrix1} {matrix2} -o {output}".format(
          matrix1 = matrix1,
          matrix2 = matrix2,
          output = output_dir+"joined_matrix.gz"))
    cmd=";".join(cbind_cmd)
    subprocess.run(cmd, shell=True)

def reorder_matrix(matrix,configfile): ##TODO
    orderedbed=os.path.join(configfile['outputDir'],"ordered.bed")
    if configfile["kmeans_clust"] or configfile["hclust"] or configfile["refIndex"]:
         matrix.group_boundaries = [0]
         ordered_regions = []
         ordered_matrix = []

         regions=zip(*matrix.regions)
         order = pd.read_csv(orderedbed, sep ='\t') #TODO what to compare to if no ordered.bed??
         from itertools import groupby, accumulate
         groups_freq = {key:len(list(group)) for key, group in groupby(order["deepTools_group"])}
         matrix.group_labels = list(groups_freq.keys())
         freq = list(accumulate(groups_freq.values()))
         matrix.group_boundaries.extend(freq)
         id_list = []
         for row in matrix.regions:
             id_list.append(row[2])
         print(len(id_list))
         match = lambda a, b: [ b.index(x) if x in b else None for x in a ]
         ii_match= match(order["name"], list(regions)[2])
         print(len(ii_match))
         ordered_regions = [ matrix.regions[i] for i in ii_match ]
         ordered_matrix = matrix.matrix[ii_match, :]

         matrix.regions = ordered_regions
         matrix.matrix = ordered_matrix
    else: #TODO
        #update group_label
        matrix.group_labels = ["genes"]

    outfile = os.path.join(configfile['outputDir'], "UpdatedclosestGene.matrix.gz")
    matrix.save_matrix(outfile)
    return outfile


def computefinalmatrix(regions, bigwigs, configfile):
    matrix_output=os.path.join(configfile['outputDir'], configfile['mode']+"_allsamples.matrix.gz")
    compute_matrix(configfile['mode'], bigwigs, regions, matrix_output, configfile)
    if configfile['extramatrix']:
       print("load matrix")
       matrix2,params2 = read_matrix_file(configfile['extramatrix'])
       print("reorder")
       additional_matrix = reorder_matrix (matrix2,configfile)
       __cbind_matrix(matrix_output,additional_matrix, configfile['outputDir'])
