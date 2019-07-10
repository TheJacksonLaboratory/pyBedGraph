from distutils.core import setup
from setuptools.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension(name="include_missing_bp",
              sources=["include_missing_bp.pyx"]),
    Extension(name="ignore_missing_bp",
              sources=["ignore_missing_bp.pyx"]),
    Extension(name="util",
              sources=["util.pyx"])
]

setup(
    ext_modules=cythonize(extensions, language_level="3")
)
