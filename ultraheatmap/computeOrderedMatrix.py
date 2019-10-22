#!/usr/bin/env python

import os
import sys
import argparse
import yaml
# import necessary modules
import ultraheatmap.computeMatrix as cm

sys.path.append(os.path.join(os.path.dirname(os.path.dirname((os.path.abspath(
                __file__))))))

# Code directory:
configDir = os.path.dirname(os.path.realpath(__file__))
# Parse entered arguments


def parse_args(defaults={}):

    """
    Parse the arguments from the command line
    """
    parser = argparse.ArgumentParser(description="The program sorts/clusters "
                                     "regions considering the reference "
                                     "samples ( given by --groupUsingSamples) "
                                     "and makes a matrix over all the samples "
                                     "using the sorted/clustered regions.",
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    # Required arguments
    required.add_argument("-S",
                          "--scoreFileName",
                          dest="bigwigs",
                          nargs='+',
                          help="bigwig files, the ordered matrix is computed"
                          "from.",
                          required=True)
    required.add_argument("-R",
                          "--regionsFileName",
                          dest="regionOfInterest",
                          nargs='+',
                          help="BED files definig the genomic regions of the "
                          "matrix.Multiple files can be provided, but the per "
                          "group information will be lost due to the clustering",
                          required=True)
    required.add_argument("-o",
                          "--outFileName",
                          dest="matrixOutput",
                          help="Matrix clustered by the given reference samples",
                          required=True)

    # optional arguments
    optional.add_argument("-g",
                          "--groupUsingSamples",
                          dest="refIndices",
                          nargs='+',
                          help="sample indices (order of the bigwig files "
                          "given via -S).It is 1-based and is used to define "
                          "the reference samples. The reference samples will "
                          "be used for sorting/clustering the regions (given "
                          "bed files), before all samples will be used to "
                          "generate the output matrix. Several indices can be "
                          "added while separated by space from each other. "
                          "Default is None and will take all the "
                          "samples into account to sort/cluster the regions.",
                          type=int,
                          default=None)

    optional.add_argument("-p",
                          "--numberOfProcessors",
                          dest="numberOfProcessors",
                          help="From deepTools doc: Number of processors to "
                          "use. Type \"max/2\" to use half the maximum number "
                          "of processors or \"max\" to use all available "
                          "processors.",
                          nargs='+',
                          type=int,
                          default=[1, 1])

    optional.add_argument("--outFileSortedRegions",
                          dest="outFileSortedRegions",
                          help="From deepTools doc: File name in which the "
                          "regions are saved after skiping zeros or min/max "
                          "threshold values. The order of the regions in the "
                          "file follows the sorting order selected. This is "
                          "useful, for example, to generate other heatmaps "
                          "keeping the sorting of the first heatmap.",
                          default=None)

    optional.add_argument("--outputReferenceMatrix",
                          dest="outputReferenceMatrix",
                          help="Matrix on the reference sampels only before "
                          "clustering",
                          default=None)

    optional.add_argument("--kmeans",
                          dest="kmeans",
                          metavar="INT",
                          type=int,
                          help="number of clusters in k-means clustering",
                          default=None)

    optional.add_argument("--hclust",
                          dest="hclust",
                          metavar="INT",
                          type=int,
                          help="Number of clusters to compute using hierarchical"
                          "clustering as defined by deepTools plotHeatmap",
                          default=None)

    optional.add_argument("-b", "--upstream",
                          "--beforeRegionStartLength",
                          dest="beforeRegionStartLength",
                          help="From deepTools doc: Distance upstream of the "
                          "start site of the regions defined in the region file."
                          " If the regions are genes, this would be the distance"
                          " upstream of the transcription start site.",
                          nargs='+',
                          type=int,
                          default=[0, 0])

    optional.add_argument("-a", "--downstream",
                          "--afterRegionStartLength",
                          dest="afterRegionStartLength",
                          help="From deepTools doc: Distance downstream of the "
                          "end site of the given regions. If the regions are "
                          "genes, this would be the distance downstream of the "
                          "transcription end site. ",
                          type=int,
                          nargs='+',
                          default=[0, 0])

    optional.add_argument("-op",
                          "--plotOutput",
                          dest="plotOutput",
                          help="File name to save the intermediate heatmap. "
                          "The file ending will be used to determine the format "
                          "of the image . Available formats are: \"png\", "
                          "\"eps\", \"pdf\" and \"svg\" (From deeptools doc)",
                          default=None)

    optional.add_argument("--config",
                          dest="userconfig",
                          help="Added to the default configuration, overwrites "
                          "if necessary.",
                          default=None)
    return parser


def merge_dictionaries(a, b):
    """
    merging two dictionaries
    """
    merged = {**a, **b}
    return merged


def add_diff(a, b):
    """
    Add the difference between the default config file and argumnets to the"
    "configfile. Order matters!
    b is defaultconfig
    a is vars(args)
    """
    for key in b.keys():
        if key in a:
            b[key] = a[key]


def check_config_file_items(configfile):
    params_list = ["regionBodyLength", "startLabel", "endLabel",
                   "unscaled5prime", "unscaled3prime", "referencePoint",
                   "nanAfterEnd", "beforeRegionStartLength",
                   "afterRegionStartLength", "binSize", "sortRegions",
                   "sortUsing", "averageTypeBins", "missingDataAsZero",
                   "skipZeros", "minThreshold", "maxThreshold",
                   "blackListFileName", "smartLabels", "verbose", "scale",
                   "numberOfProcessors", "metagene", "transcriptID", "exonID",
                   "transcript_id_designator"]
    for item in params_list:
        assert len(configfile[item]) <= 2, "length of "+item+" should be <= 2!"
        if len(configfile[item]) == 1:
            configfile[item].append(configfile[item][0])


def main():
    """
    Compute a matrix from an ordered region file
    """
    # First part of the code is applied for all cases
    defaultconfigfile = {}
    # 1. Read the default config file
    with open(os.path.join(configDir, "configs", "computeOrderedMatrix.yaml"),
              'r') as stream:
        defaultconfigfile = yaml.load(stream)

    # 2. Parse the arguments
    parser = parse_args(defaultconfigfile)
    args = parser.parse_args()

    # 2.1 modify config file if user config is provided
    configfile = defaultconfigfile
    if args.userconfig:
        with open(os.path.join(args.userconfig), 'r') as stream:
            userconfigfile = yaml.load(stream)
            configfile = merge_dictionaries(configfile, userconfigfile)

    configfile = merge_dictionaries(configfile, vars(args))

    # 3. check the length of computeMatrix arguments to be equal to 2
    check_config_file_items(configfile)

    # 4. Generate an ordered region, using references only
    if configfile["outFileSortedRegions"] is None:
        path_name = os.path.dirname(os.path.abspath(args.matrixOutput))
        configfile["outFileSortedRegions"] = path_name+'/orderedBedFile.bed'

    if configfile["refIndices"]is None:
        configfile["refIndices"] = list(range(1, len(configfile["bigwigs"])+1))
    cm.sortbyreference(configfile["regionOfInterest"], configfile["bigwigs"],
                       configfile["refIndices"], configfile)
    assert(os.path.getsize(configfile["outFileSortedRegions"]) > 0)

    # 4.Build a matrix over all the samples
    hm = cm.computefinalmatrix(configfile["outFileSortedRegions"],
                               configfile["bigwigs"], configfile)

    matrix_output = os.path.join(args.matrixOutput)
    hm.save_matrix(matrix_output)
