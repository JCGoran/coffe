#!/usr/bin/env python3

from setuptools import Extension, setup
from Cython.Build import cythonize

setup(
    name="Coffe",
    version="3.0",
    description="The COrrelation Function Full-sky Estimator code",
    url="https://github.com/JCGoran/coffe",
    author="Goran Jelic-Cizmek",
    author_email="goran.jelic-cizmek@unige.ch",
    install_requires=["numpy>=1.19.5"],
    packages=["coffe"],
    ext_modules=cythonize([Extension("coffe.coffe", ["coffe/*.pyx"])]),
)
