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

             git clone https://github.com/maxplanck-ie/ultraheatmap.git

Afterwards, create a new conda environment with all the Prerequisites by running the following command line:

              conda env create -f requirements.yaml


Then activate the environment:

              conda activate ultraheatmap

To install the program in this environment:

              python setup.py install

Alternatively, `pip` or `conda` can be used to install the package. We highly recommend you to create a new conda environment prior to the installation and install it after activating this environment.

from the ultraheatmap directory. Now, you already have the program installed and can access each of the modules by calling them. Try

            $  ultraheatmap -h ,

            $  computeOrderedMatrix -h

or


            $  AddFeatureToMatrix -h

To terminate the environment run:

              conda deactivate
