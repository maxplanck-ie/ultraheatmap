from distutils.core import setup

setup(
    name='ultraHeatmap',
    packages=['computeOrderedMatrix', 'BuildDeeptoolsLikeMatrix'],
    long_description=open('README.md').read(),
    install_requires=[ #TODO
        "deeptools=3.0.1",
        "pybedtools=0.7.10",
        "gffutils=0.9",
        "numpy=1.13.3",
        "pandas=0.22.0",
        "yaml=0.1.7"
    ],
    cmdclass={'sdist': sdist, 'install': install}
)
