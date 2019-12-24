import setuptools
from distutils.extension import Extension

USE_CYTHON = False
try:
    from Cython.Distutils import build_ext
    USE_CYTHON = True
except ImportError:
    pass

ext = '.pyx' if USE_CYTHON else '.c'

extensions = [
    Extension('pyBedGraph.ignore_missing_bp', ['pyBedGraph/ignore_missing_bp' + ext]),
    Extension('pyBedGraph.include_missing_bp', ['pyBedGraph/include_missing_bp' + ext]),
    Extension('pyBedGraph.util', ['pyBedGraph/util' + ext]),
]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions, language_level=3)

NAME = 'pyBedGraph'
VERSION = '0.5.41'

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

    install_requires=['numpy>=1.16.4'],

    python_requires='>=3.6',

    ext_modules=extensions,

    data_files=[("", ["LICENSE"])],

    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],

)
