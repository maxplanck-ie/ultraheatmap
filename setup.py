# -*- coding: utf-8 -*-

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
    name='ultraheatmap',
    packages=find_packages(),
    scripts=['bin/computeOrderedMatrix', 'bin/AddFeatureToMatrix',
             'bin/ultraheatmap'],
    license="MIT",
    long_description=open('README.md').read(),
    include_package_data=True,
    cmdclass={'sdist': sdist, 'install': install}
)
