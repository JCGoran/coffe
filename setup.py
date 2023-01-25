#!/usr/bin/env python3

import os
import subprocess

from Cython.Build import cythonize
from setuptools import Extension, setup


def get_include_dirs():
    """
    Returns a list of all of the alternative include directories (virtual env,
    conda env, etc.)
    """
    return [
        os.path.join(os.environ.get(_), "include/")
        for _ in ["VIRTUAL_ENV", "CONDA_PREFIX"]
        if os.environ.get(_)
    ]


def get_library_dirs():
    """
    Returns a list of all of the alternative library directories (virtual env,
    conda env, etc.)
    """
    return [
        os.path.join(os.environ.get(_), "lib/")
        for _ in ["VIRTUAL_ENV", "CONDA_PREFIX"]
        if os.environ.get(_)
    ]


extra_compile_args = [
    "-fopenmp",
    "-Ofast",
    "-DHAVE_CLASS",
    "-DHAVE_CUBA",
    "-DCOFFE_CYTHON",
]


commit = subprocess.run(
    ["git", "rev-parse", "HEAD"],
    capture_output=True,
    encoding="utf-8",
).stdout.strip()

setup(
    name="coffe",
    version="3.0.0",
    description=commit,
    url="https://github.com/JCGoran/coffe",
    author="Goran Jelic-Cizmek",
    author_email="goran.jelic-cizmek@unige.ch",
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
                ]
                + get_include_dirs(),
                libraries=[
                    "m",
                    "gsl",
                    "gslcblas",
                    "fftw3",
                    "cuba",
                    "class",
                ],
                library_dirs=get_library_dirs(),
                extra_compile_args=extra_compile_args,
                extra_link_args=["-fopenmp"],
            ),
        ]
    ),
)
