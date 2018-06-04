import subprocess
import shlex
import sys
import os

sys.path.append(os.path.join(os.path.dirname(os.path.dirname((os.path.realpath(__file__))))))
from bashwrapper.bashwrapper import Bash
from deeptoolsapi.plotHeatMap import plot_heatmap
from deeptools.heatmapper import heatmapper as dh #TODO
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


def __cbind_matrix(matrix1, matrix2):
    cbind_cmd = [deeptools_module_load]
    cbind_cmd.append("computeMatrixOperations cbind -m {matrix1} {matrix2} -o {output}".format(
          matrix1 = matrix1,
          matrix2 = matrix2,
          output = "joined_matrix.gz"))
    cmd=";".join(cbind_cmd)
    subprocess.run(cmd, shell=True)


def computefinalmatrix(regions, bigwigs, configfile):
    matrix_output=os.path.join(configfile['outputDir'], configfile['mode']+"_allsamples.matrix")
    compute_matrix(configfile['mode'], bigwigs, regions, matrix_output, configfile)
    if configfile['extramatrix']:
       __cbind_matrix(matrix_output,configfile['extramatrix'])


