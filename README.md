# ultraheatmap

ultraheatmaps facilitates the production of [deepTools](https://github.com/deeptools/deepTools)
heatmaps. The heatmaps typically show signal at genomic regions, which can be
appended by orthogonal data, like associated gene expression. ultraheatmap
facilitates adding orthogonal data to a deepTools matrix and allows to cluster a
genomic heatmap by selected samples in just one single command-line call.


## Getting Started

These instructions will get you a copy of ultraheatmaps up and running on your local machine.

### Prerequisites

The program is written and tested in python3 and requires the following packages
```
import argparse
import csv
import deeptoolsapi
import gffutils
import numpy
import os
import pandas
import pybedtools
import subprocess
import sys
import unittest
import yaml
```

For installation of the following packages, please check the project documentation

* [deeptoolsapi](https://github.com/deeptools/deepTools)
* [pybedtools](https://daler.github.io/pybedtools/)  
* [gffutils](https://pythonhosted.org/gffutils/installation.html)


### Installation

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags).

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
