from distutils.core import setup
from setuptools.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension(name="complete_stats",
              sources=["complete_stats.pyx"]),
    Extension(name="ignore_blank_stats",
              sources=["ignore_blank_stats.pyx"])
]

setup(
    ext_modules=cythonize(extensions, language_level="3")
)
