import setuptools
from distutils.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension('pyBedGraph.ignore_missing_bp', ['pyBedGraph/ignore_missing_bp.pyx']),
    Extension('pyBedGraph.include_missing_bp', ['pyBedGraph/include_missing_bp.pyx']),
    Extension('pyBedGraph.util', ['pyBedGraph/util.pyx']),
]

NAME = 'pyBedGraph'
VERSION = '0.5.32'

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

    install_requires=['numpy>=1.16.4', 'pyBigWig>=0.3.16', 'cython>=0.29.12'],

    ext_modules=cythonize(extensions, language_level=3),

    data_files=[("", ["LICENSE"])],

    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],

)
