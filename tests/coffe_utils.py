"""
Various helper utilities that are outside of the scope of the core functionality of COFFE.
"""

from __future__ import annotations

from typing import Callable, Optional, Union

import numpy as np
import pandas as pd
from coffe import AverageCovariance, Covariance
from scipy.linalg import block_diag


def covariance_matrix(
    cov: list[Covariance],
    l: Optional[list[int]] = None,
    z_mean: Optional[list[float]] = None,
    deltaz: Optional[list[float]] = None,
    rmin: Optional[Union[Callable, float, int]] = None,
    rmax: Optional[Union[Callable, float, int]] = None,
    rstep: Optional[Union[float, int]] = None,
):
    """
    Converts an array of covariances into a numpy matrix for easy matrix
    multiplication.

    Parameters
    ----------
    cov : array_like of covariances
        the covariances, in the same format as output of
        `coffe.Coffe.compute_covariance_bulk`

    l : Optional[list[int]], default = None
        the list of values of l

    z_mean : Optional[list[float]], default = None
        the list of values of mean redshifts

    deltaz : Optional[list[float]], default = None
        the list of redshift bin widths. Only useful if one of `rmin`, `rmax` is a function.

    rmin : Optional[Union[Callable, float, int]], default = None
        If an int or a float, the smallest separation to take from the covariance matrix.
        If a function, must take exactly two mandatory arguments, which should
        be the mean redshift and the size of the bin.
        An example can be `lambda z_mean, deltaz: 1 / np.sqrt(1 + z)`

    rmax : Optional[Union[Callable, float, int]], default = None
        If an int or a float, the largest separation to take from the covariance matrix.
        If a function, must take exactly two mandatory arguments, which should
        be the mean redshift and the size of the bin.
        An example is the function `maximum_separation` from the `Coffe` class.

    rstep : Optional[float] = None
        if set, only takes separations from the covariance matrix that are
        multiples of `rstep`.

    Examples
    --------
    >>> covariance_matrix(Coffe(has_density=True).compute_covariance_bulk())
    """

    def convert_array_to_matrix(arr):
        size = round(np.sqrt(len(arr)))
        return np.reshape(arr, (size, size))

    if not all(isinstance(_, Covariance) for _ in cov):
        raise ValueError(f"{cov} is not an array of covariances")

    df = pd.DataFrame([_.to_dict() for _ in cov])

    # N.B. this may require fixing due to floating point math
    if z_mean is not None:
        df = df.loc[df.z.isin(z_mean)]
    else:
        z_mean = df.z.unique()

    if l is not None:
        df = df.loc[(df.l1.isin(l)) & (df.l2.isin(l))]

    if rmin is not None:
        if isinstance(rmin, (int, float)):
            df = df.loc[(df.r1 >= rmin) & (df.r2 >= rmin)]
        else:
            func_rmin = rmin
            df_rmin = pd.DataFrame()
            for z, dz in zip(z_mean, deltaz):
                df_rmin = df_rmin.append(
                    df.loc[
                        (df.r1 >= func_rmin(z, dz))
                        & (df.r2 >= func_rmin(z, dz))
                        & (df.z == z)
                    ]
                )
            df = df_rmin

    if rmax is not None:
        if isinstance(rmax, (int, float)):
            df = df.loc[(df.r1 <= rmax) & (df.r2 <= rmax)]
        else:
            func_rmax = rmax
            df_rmax = pd.DataFrame()
            for z, dz in zip(z_mean, deltaz):
                df_rmax = df_rmax.append(
                    df.loc[
                        (df.r1 <= func_rmax(z, dz))
                        & (df.r2 <= func_rmax(z, dz))
                        & (df.z == z)
                    ]
                )
            df = df_rmax

    if rstep is not None:
        df = df[np.isclose(df.r1 % rstep, 0)]

    return block_diag(
        *[
            convert_array_to_matrix(df.loc[df.z == z].value.to_numpy(dtype=float))
            for z in z_mean
        ]
    )


def average_covariance_matrix(
    cov: list[AverageCovariance],
    l: Optional[list[int]] = None,
    z_min: Optional[list[float]] = None,
    z_max: Optional[list[float]] = None,
    rmin: Optional[Union[Callable, float, int]] = None,
    rmax: Optional[Union[Callable, float, int]] = None,
    rstep: Optional[Union[float, int]] = None,
):
    """
    Converts an array of covariances into a numpy matrix for easy matrix
    multiplication.

    Parameters
    ----------
    cov : array_like of covariances
        the covariances, in the same format as output of
        `coffe.Coffe.compute_average_covariance_bulk`

    l : Optional[list[int]], default = None
        the list of values of l

    z_min : Optional[list[float]], default = None
        the list of values of min redshifts

    z_max : Optional[list[float]], default = None
        the list of values of max redshifts

    rmin : Optional[Union[Callable, float, int]], default = None
        If an int or a float, the smallest separation to take from the covariance matrix.
        If a function, must take exactly two mandatory arguments, which should
        be the mean redshift and the size of the bin.
        An example can be `lambda z_mean, deltaz: 1 / np.sqrt(1 + z)`

    rmax : Optional[Union[Callable, float, int]], default = None
        If an int or a float, the largest separation to take from the covariance matrix.
        If a function, must take exactly two mandatory arguments, which should
        be the mean redshift and the size of the bin.
        An example is the function `maximum_separation` from the `Coffe` class.

    rstep : Optional[float] = None
        if set, only takes separations from the covariance matrix that are
        multiples of `rstep`.

    Examples
    --------
    >>> average_covariance_matrix(Coffe(has_density=True).compute_average_covariance_bulk())
    """

    def convert_array_to_matrix(arr):
        size = round(np.sqrt(len(arr)))
        return np.reshape(arr, (size, size))

    if not all(isinstance(_, AverageCovariance) for _ in cov):
        raise ValueError(f"{cov} is not an array of average covariances")

    df = pd.DataFrame([_.to_dict() for _ in cov])

    # N.B. this may require fixing due to floating point math
    if z_min is not None:
        df = df.loc[df.z_min.isin(z_min)]
    else:
        z_min = df.z_min.unique()

    if z_max is not None:
        df = df.loc[df.z_max.isin(z_max)]
    else:
        z_max = df.z_max.unique()

    if l is not None:
        df = df.loc[(df.l1.isin(l)) & (df.l2.isin(l))]

    return block_diag(
        *[
            convert_array_to_matrix(df.loc[df.z_min == z].value.to_numpy(dtype=float))
            for z in z_min
        ]
    )
