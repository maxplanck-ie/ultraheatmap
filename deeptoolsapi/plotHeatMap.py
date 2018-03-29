import subprocess
import os

###Load deeptools
deeptools_module_load="module load deeptools"

#######plotHeatMap
def plot_heatmap(mode, matrix, configfile, ordered):
   """
   """
   heatmap_cmd =[deeptools_module_load]
   heatmap_plot = os.path.join(configfile['outputDir'],mode+"_heatmap.png")
   kmeans=""
   hclust=""
   if configfile['kmeans_clust'] != 0:
      kmeans=' --kmeans '+ str(configfile['kmeans_clust'])

   if configfile["hclust"] != 0:
      hclust=' --hclust '+str(configfile['hclust'])

   if ordered == True:
          heatmap_plot = os.path.join(configfile['outputDir'],mode+"_refonly_heatmap.png")
          orderedbed = os.path.join(configfile['outputDir'],"ordered.bed")
          heatmap_cmd.append("plotHeatmap -m "+matrix+" --outFileName "+ heatmap_plot + " --outFileSortedRegions "+orderedbed+kmeans+hclust)
   else:
          print("heat_map on all samples")
          heatmap_cmd.append("plotHeatmap -m "+matrix+" --outFileName "+ heatmap_plot + " --sortRegions no "+kmeans+hclust)

   cmd=";".join(heatmap_cmd)
   print(cmd)
   subprocess.run(cmd, shell=True)

