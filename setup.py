#from distutils.core import setup
from setuptools import setup, Extension, find_packages
from setuptools.command.sdist import sdist as _sdist
from setuptools.command.install import install as _install

class sdist(_sdist):

    def run(self):
        return _sdist.run(self)


class install(_install):

    def run(self):
        _install.run(self)
        return

setup(
    name='ultraHeatmap',
    packages=find_packages(),
    scripts=['bin/computeOrderedMatrix', 'bin/BuildADeeptoolsLikeMatrix'],
    long_description=open('README.md').read(),
    # install_requires=[ #TODO
    #     "deeptools=3.0.1",
    #     "pybedtools=0.7.10",
    #     "gffutils=0.9",
    #     "numpy=1.13.3",
    #     "pandas=0.22.0",
    #     "yaml=0.1.7"
    # ],
     cmdclass={'sdist': sdist, 'install': install}
)
