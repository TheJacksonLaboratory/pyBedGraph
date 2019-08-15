import setuptools
from distutils.core import setup
from distutils.extension import Extension

extensions = [
    Extension('pyBedGraph.ignore_missing_bp', ['pyBedGraph/ignore_missing_bp.c']),
    Extension('pyBedGraph.include_missing_bp', ['pyBedGraph/include_missing_bp.c']),
    Extension('pyBedGraph.util', ['pyBedGraph/util.c']),
]

cmdclass = {}

NAME = 'pyBedGraph'
VERSION = '0.5.22'

setuptools.setup(

    name=NAME,

    version=VERSION,

    author="Henry Zhang",

    author_email="henrybzhang.99@gmail.com",

    description="A package for fast operations on 1-dimensional genomic signal tracks",

    long_description=open('README.md').read(),

    long_description_content_type="text/markdown",

    url="https://github.com/TheJacksonLaboratory/pyBedGraph",

    packages=setuptools.find_packages(),

    install_requires=['numpy', 'pyBigWig'],

    ext_modules=extensions,

    cmdclass=cmdclass,

    classifiers=[

        "Programming Language :: Python :: 3",

        "License :: OSI Approved :: MIT License",

        "Operating System :: OS Independent",

    ],

)
