# generating hybrid matrices and plot them by plotHeatmap

#!/usr/bin/env python

import os
import subprocess as sp
import glob
import sys
import argparse
import yaml

def parse_args():
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    # required argumnets:
    required.add_argument("--config",
                          "-c",
                          dest="config",
                          help="config to use to plot heatmap (yaml file)",
                          required=True)

    # optional arguments:
    optional.add_argument("--heatmapOnly",
                           "-ho",
                           default=False,
                           action='store_true',
                           dest="heatmapOnly",
                           help="only plot heatmap for given matrices in yaml file")



    return parser


def main():
    """
    compute_matrix > ultraheatmap (optional) > heatmap (or profile only if no ultraheatmap)
    """
    parser = parse_args()
    args = parser.parse_args()
    with open(os.path.join(args.config), 'r') as stream:
        config = yaml.safe_load(stream)


    # transcriptomics data
    try:
        path_2_rna = config["path_2_rna"]
    except:
        path_2_rna = ""

    # regions to plot on
    regions_to_plot = config["regions_to_plot"]

    # output params
    output_path = config["output_path"]

    matrix_name = os.path.join(output_path, config["matrix_name"])
    scale_region = config["scale_region"]

    heatmap_name = os.path.join(output_path, config["heatmap_name"])

    bws = " "
    names = " "

    for k, v in config["bws"].items():
        bws += v+" "
        names += k+" "

    deg_list =[]
    deg_names = " "
    if config["ultraheatmap_matrix_name"] != "":
        for i, deg in enumerate( config["deg_list"]):
            deg_list.append(os.path.join(config["path_2_rna"], deg))
            # deg_list = [os.path.join(path_2_rna, "nej_vs_wt_shrunk.tsv")]
            deg_names += config["deg_names"][i]+" "
            ultraheatmap_matrix_name = os.path.join(output_path, config["ultraheatmap_matrix_name"])


    if not args.heatmapOnly:
        if scale_region:
            cmd = "computeMatrix scale-regions "
            cmd += " -S "+bws
            cmd += " -R "+regions_to_plot
            cmd += " -p 20 "
            cmd += "-bs "+str(config["binsize"])
            cmd += " --skipZeros --regionBodyLength "+str(config["regionBodyLength"])
            cmd += " -a "+str(config["after_region"])
            cmd += " -b "+str(config["before_region"])
            cmd += " --missingDataAsZero --skipZeros --samplesLabel "+names
            cmd += " -o "+ matrix_name
            if config['averageTypeBins']:
                cmd += ' --averageTypeBins '+config['averageTypeBins']
            print(cmd)
            sp.check_output(cmd, shell = True)
        else:
            cmd = "computeMatrix reference-point --referencePoint "+config["refpoint"]
            cmd += " -S "+bws
            cmd += " -R "+regions_to_plot
            cmd += " -p 20 --skipZeros "
            cmd += " -a "+str(config["after_region"])
            cmd += " -b "+str(config["before_region"])
            cmd += " --missingDataAsZero --skipZeros --samplesLabel "+names
            cmd += " -o "+ matrix_name
            if config['averageTypeBins']:
                cmd += ' --averageTypeBins '+config['averageTypeBins']
            sp.check_output(cmd, shell = True)
        if config["ultraheatmap_matrix_name"] != "":
            w_mapping = config["w_mapping"]
            if w_mapping:
                cmd = "addFeatureToMatrix "
                cmd += "-m "+ matrix_name
                cmd += " -o "+ultraheatmap_matrix_name
                cmd += " -t "+" ".join(deg_list)
                cmd += " --featureNames  "+config['featureNames']
                cmd += " --referencePoint TSS "
                cmd += " --genomeGtf "+config["gtf"]
                sp.check_output(cmd, shell = True)
            else:
                cmd = "addFeatureToMatrix "
                cmd += "-m "+ matrix_name
                cmd += " -o "+ultraheatmap_matrix_name
                cmd += " -t "+" ".join(deg_list)
                if config['featureNames']:
                    cmd += " --featureNames "+config['featureNames']
                if config['featureIdColumn']:
                    cmd += " --featureIdColumn "+config['featureIdColumn']
                print(cmd)
                sp.check_output(cmd, shell = True)


    # plot heatmap
    all_samples_names = names
    matrix_to_plot = matrix_name
    if config["ultraheatmap_matrix_name"] != "":
        all_samples_names+= deg_names
        matrix_to_plot = ultraheatmap_matrix_name
    cmd = "plotHeatmap -m "
    cmd += matrix_to_plot
    if config["zmin"]:
        cmd += " --zMin  "+config["zmin"]
    if config["zmax"]:
        cmd+= " --zMax  "+config["zmax"]
    if config["colors"]:
        cmd += " --colorMap "+config["colors"]
    cmd += " -o "+heatmap_name
    cmd += " --samplesLabel "+all_samples_names
    cmd += " --whatToShow 'heatmap and colorbar' "
    if sort_by_sample != "":
        cmd += " --sortUsingSamples "+config["sort_by_sample"]
    if config["sortUsing"]:
        cmd += " --sortUsing "+config["sortUsing"]
    if config["sorting_direction"]
        "--sortRegions "+config["sorting_direction"]
    if config["sorted_regions"] != "":
        sorted_regions = os.path.join(output_path, config["sorted_regions"])
        cmd += " --outFileSortedRegions "+sorted_regions
    if config["sorted_matrix"] != "":
        sorted_matrix = os.path.join(output_path, config["sorted_matrix"])
        cmd += "  --outFileNameMatrix "+sorted_matrix
    if config["kmeans"]:
        cmd += " --kmeans "+config["kmeans"]
        if config["clusterUsingSamples"]:
            cmd += " --clusterUsingSamples "+config["clusterUsingSamples"]
    else:
        if config["reg_labels"]:
            cmd += " --regionsLabel "+config["reg_labels"]
    print(cmd)
    sp.check_output(cmd, shell = True)
if __name__ == "__main__":
   main()
