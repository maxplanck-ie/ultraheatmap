# ultraheatmap

ultraheatmap facilitates the production of [deepTools](https://github.com/deeptools/deepTools)
heatmaps. The heatmaps typically show signal at genomic regions, which can be
appended by orthogonal data, like associated gene expression. ultraheatmap
facilitates adding orthogonal data to a deepTools matrix and allows to cluster a
genomic heatmap by selected samples in just one single command-line call.


## Getting Started

These instructions will get you a copy of ultraheatmap up and running on your local machine.

### Prerequisites

The Prerequisites can be found in requirements.yaml

### Installation

First, get the source code:

             $  git clone https://github.com/maxplanck-ie/ultraheatmap.git

Afterwards, create a new conda environment with all the Prerequisites by running the following command line:


              $ conda env create -f requirements.yaml


Then activate the environment:

              $ conda activate ultraheatmap

To install the program in this environment:

              $ python setup.py install
from the ultraheatmap directory.

Alternatively, `pip` or `conda` can be used to install the package. We highly recommend you to create a new conda environment prior to the installation and install it after activating this environment. This can be done as follows:

              $ conda create -n ultraheatmap python=3.6

              $ conda activate ultraheatmap

              $ conda install -c bioconda -c conda-forge ultraheatmap


Now, you already have the program installed and can access each of the modules by calling them. Try

              $  ultraheatmap -h ,

              $  computeOrderedMatrix -h

or

              $  AddFeatureToMatrix -h

To terminate the environment run:

              $ conda deactivate

usage example:

  $ computeOrderedMatrix -S bigwig1.bw bigwig2.bw -R 2peaks.bed -o final_matrix.gz -p 20 -a 100 --outputReferenceMatrix \     intermediate_matrix.gz -op intermediate_matrix_heatmap.png -g 1 --kmeans 2

  the above command line produces a `deeptools` matrix on both given bw while the regions are the clusters obtained from the given bed files after using kmeans clustering algorithm with 2 clusters (`--kmeans 2`) based on the signal of first bigwig file (`-g 1`).
