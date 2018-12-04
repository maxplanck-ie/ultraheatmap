import os
import sys
import numpy as np

#######plotHeatMap
def __plot_heatmap(hm,indexList, configfile):
   """
   """
   if hm.parameters['min threshold'] is not None or hm.parameters['max threshold'] is not None:
        hm.filterHeatmapValues(hm.parameters['min threshold'], hm.parameters['max threshold'])

   if configfile["sortRegions"] == 'keep':
      configfile["sortRegions"] = 'no'  # These are the same thing ##XXX what does that mean???

   if configfile["kmeans"] is not None:
        hm.matrix.hmcluster(configfile["kmeans"], method='kmeans')
   else:
        if configfile["hclust"] is not None:
            print("Performing hierarchical clustering."
                  "Please note that it might be very slow for large datasets.\n")
            hm.matrix.hmcluster(configfile["hclust"], method='hierarchical')

   group_len_ratio = np.diff(hm.matrix.group_boundaries) / len(hm.matrix.regions)
   if np.any(group_len_ratio < 5.0 / 1000):
        problem = np.flatnonzero(group_len_ratio < 5.0 / 1000)
        sys.stderr.write("WARNING: Group '{}' is too small for plotting, you might want to remove it. "
                         "There will likely be an error message from matplotlib regarding this "
                         "below.\n".format(hm.matrix.group_labels[problem[0]]))


   if configfile["sortRegions"] != 'no':
        hm.matrix.sort_groups(sort_using=configfile["sortUsing"],sort_method=configfile["sortRegions"],sample_list=indexList)

   assert(configfile["outFileSortedRegions"])
   hm.save_BED(open(configfile["outFileSortedRegions"], "w"))
