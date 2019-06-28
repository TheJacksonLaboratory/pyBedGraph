import setuptools

NAME = 'pyBedGraph'
VERSION = '0.3.1'

with open("README.md", "r") as fh:

    long_description = fh.read()

setuptools.setup(

    name=NAME,

    version=VERSION,

    author="Henry Zhang",

    author_email="henrybzhang.99@gmail.com",

    description="An alternative to pyBigWig for bedgraph files",

    long_description=long_description,

    long_description_content_type="text/markdown",

    url="https://github.com/c0ver/pyBedGraph",

    packages=setuptools.find_packages(),

    classifiers=[

         "Programming Language :: Python :: 3",

         "License :: OSI Approved :: MIT License",

         "Operating System :: OS Independent",

    ],

 )
