#!/usr/bin/env python3

import os

from Cython.Build import cythonize
from setuptools import Extension, setup

extra_compile_args = [
    "-fopenmp",
    "-Ofast",
    "-DHAVE_CLASS",
    "-DHAVE_CUBA",
    "-DCOFFE_CYTHON",
]

setup(
    name="Coffe",
    version="3.0",
    description="The COrrelation Function Full-sky Estimator code",
    url="https://github.com/JCGoran/coffe",
    author="Goran Jelic-Cizmek",
    author_email="goran.jelic-cizmek@unige.ch",
    install_requires=["numpy>=1.19.5"],
    packages=["coffe"],
    ext_modules=cythonize(
        [
            Extension(
                "coffe.coffe",
                sources=[
                    "coffe/*.pyx",
                    "src/errors.c",
                    "src/common.c",
                    "src/parser.c",
                    "src/background.c",
                    "src/twofast.c",
                    "src/integrals.c",
                    "src/signal.c",
                    "src/functions.c",
                    "src/corrfunc.c",
                    "src/multipoles.c",
                    "src/utils.c",
                    "src/twobessel.c",
                    "src/covariance.c",
                ],
                include_dirs=[
                    "src/",
                    "./",
                ],
                libraries=[
                    "m",
                    "gsl",
                    "gslcblas",
                    "fftw3",
                    "cuba",
                    "class",
                ],
                extra_compile_args=extra_compile_args,
                extra_link_args=["-fopenmp"],
            ),
        ]
    ),
)
