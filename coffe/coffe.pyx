# cython: binding=True
# cython: language_level=3

# TODO figure out how to use OpenMP

from libc.stdlib cimport malloc, free
from coffe cimport ccoffe
from coffe.representation import (
    Corrfunc,
    Multipoles,
    Covariance,
)

from cython.operator cimport dereference
from ctypes import CFUNCTYPE

import os
import sys
import warnings
from configparser import ConfigParser
from pathlib import Path
from typing import Any, Callable, List, Tuple, Union, Optional, Sequence

import numpy as np


class LegacyOption:
    """
    Class for handling legacy options in the config file.
    """

    def __init__(self, name, kind, replaced_by = None):
        """
        Constructor

        Parameters
        ----------
        name : str
            Name of the legacy option

        kind : str
            Type of the legacy option

        replaced_by : str, optional
            Name of the new option that replaces the legacy one
        """
        self.name = name
        self.kind = kind
        self.replaced_by = replaced_by



class CoffeConfigParser(ConfigParser):
    """
    The custom parser for the configuration file of COFFE
    """

    def getfloat_array(self, section, option, *args, **kwargs):
        """
        Converts the value of the option to a list of floats
        """
        return [float(_) for _ in self.get(section, option).strip("[]").split(",")]

    def getint_array(self, section, option, *args, **kwargs):
        """
        Converts the value of the option to a list of ints
        """
        return [int(_) for _ in self.get(section, option).strip("[]").split(",")]

    def getstr_array(self, section, option, *args, **kwargs):
        """
        Converts the value of the option to a list of strings
        """
        return [
            str(_.strip(' "')) for _ in self.get(section, option).strip("[]").split(",")
        ]



def _check_parameter(
    name : str,
    value : Any,
    kind : Union[type, Tuple[type], List[type]],
    xmin : float = None,
    xmax : float = None
):
    """
    Performs type and value checking for the parameter called `name`.
    """
    if not isinstance(value, kind):
        raise TypeError(
            f'Expected {kind}, got {type(value)}'
        )

    if xmin is not None and value < xmin:
        raise ValueError(
            f'Value {value} smaller than xmin ({xmin})'
        )

    if xmax is not None and value > xmax:
        raise ValueError(
            f'Value {value} larger than xmax ({xmax})'
        )



cdef int get_spline_size(
    ccoffe.coffe_interpolation *interp,
) except *:
    return dereference(interp.spline).size



cdef double get_spline_min(
    ccoffe.coffe_interpolation *interp,
) except *:
    return dereference(interp.spline).x[0]



cdef double get_spline_max(
    ccoffe.coffe_interpolation *interp,
) except *:
    return dereference(interp.spline).x[get_spline_size(interp) - 1]



cdef double evaluate_spline(
    ccoffe.coffe_interpolation *interp,
    const double x,
) except *:
    """
    Evaluates the spline at some value x
    """
    x_min = dereference(interp.spline).x[0]
    x_max = dereference(interp.spline).x[dereference(interp.spline).size - 1]

    _check_parameter('x', x, (int, float), x_min, x_max)

    return ccoffe.coffe_interp_spline(interp, x)



cdef int set_spline(
    ccoffe.coffe_interpolation *interp,
    x_sampling : List[float],
    y_sampling : List[float],
    const ccoffe.coffe_interp1d_type interp1d_type,
    xmin : float = 0,
    xmax : float = 15,
) except *:
    """
    Sets the value of the spline according to some x_sampling and y_sampling
    """
    if not np.allclose(x_sampling, np.sort(x_sampling)):
        raise ValueError('The input array must be sorted')
    if len(x_sampling) != len(y_sampling):
        raise ValueError(
            f'Mismatching lengths for x ({len(x_sampling)}) and y ({len(y_sampling)})'
        )

    _check_parameter('xmin', x_sampling[0], (int, float), xmin=xmin)
    _check_parameter('xmax', x_sampling[-1], (int, float), xmax=xmax)

    cdef double *x = NULL
    cdef double *y = NULL

    size = len(x_sampling)
    x = <double *>malloc(sizeof(double) * size)
    y = <double *>malloc(sizeof(double) * size)

    for i in range(size):
        x[i] = x_sampling[i]
        y[i] = y_sampling[i]

    ccoffe.coffe_init_spline(interp, x, y, size, interp1d_type)

    free(x)
    free(y)



_allowed_pk_types = {
    ccoffe.COFFE_PK_LINEAR : 'linear',
    ccoffe.COFFE_PK_LINEAR_CLASS : 'linear_class',
    ccoffe.COFFE_PK_NONLINEAR_HALOFIT : 'halofit',
    ccoffe.COFFE_PK_NONLINEAR_HMCODE : 'hmcode',
}

_allowed_pk_types_inverse = {
    value : key for key, value in _allowed_pk_types.items()
}



cdef class Coffe:
    r"""
    Class for handling the COFFE library.

    ## Examples

    First, import the `Coffe` class:
    >>> from coffe import Coffe

    ### Initialization

    To initialize the COFFE class with the default parameters, run:
    >>> cosmology = Coffe()

    You can also read a configuration file using `from_file`:
    >>> Coffe.from_file([FILE])

    There is of course an inverse, `to_file`, which can be used to save the
    current configuration to a file:
    >>> cosmology.to_file([FILE])

    ### Setting parameters

    You can change any of the parameters when creating the instance by passing
    them as keyword arguments:
    >>> cosmology = Coffe(omega_m=0.35)

    or later:
    >>> cosmology.omega_m = 0.32

    For a list of settable parameters, you can print the parameters:
    >>> cosmology.parameters

    You can also set the parameters in bulk using `set_parameters`:
    >>> cosmology.set_parameters(h=0.7, n_s=0.9)

    The only values which cannot be set as above are:
    * the galaxy, magnification, and evolution bias
    * the input power spectrum (**NOTE**: COFFE by default uses
    [CLASS](https://github.com/lesgourg/class_public/) to generate the linear
    matter power spectrum on-the-fly; this can be overridden if necessary)

    ### Setting bias and power spectrum parameters

    To set the biases, you can run:
    >>> cosmology.set_galaxy_bias1([REDSHIFTS], [VALUES])

    where `[REDSHIFTS]` is an array of redshifts (in increasing order) and
    `[VALUES]` are the corresponding values of the bias.
    There are several functions available for accomplishing this, all requiring
    the same kind of input as the above:
    * `set_galaxy_bias1`
    * `set_galaxy_bias2`
    * `set_galaxy_bias3`
    * `set_galaxy_bias4`
    * `set_magnification_bias1`
    * `set_magnification_bias2`
    * `set_evolution_bias1`
    * `set_evolution_bias2`

    To obtain the values of the biases, you can run:
    >>> cosmology.galaxy_bias1([REDSHIFT])

    where `[REDSHIFT]` is the redshift (must be a Python `int` or a `float`) at
    which we want to evaluate the given bias.

    COFFE can also evaluate the linear power spectrum (internally uses CLASS to
    generate it) at a given wavenumber and redshift using `power_spectrum`:
    >>> cosmology.power_spectrum([WAVENUMBER], [REDSHIFT])

    The power spectrum can also be set using `set_power_spectrum_linear`:
    >>> cosmology.set_power_spectrum_linear([WAVENUMBERS], [VALUES])

    where `[VALUES]` denote the values of the power spectrum at redshift zero.

    ### Background quantities

    There are several background (i.e. only redshift-dependent) quantities
    which can be computed with COFFE:
    * `scale_factor`: evaluates the scale factor $a(z)$
    * `hubble_rate`: evaluates the Hubble rate $H(z)$
    * `hubble_rate_conformal`: evaluates the conformal Hubble rate
    $\mathcal{H}(z)$
    * `hubble_rate_conformal_derivative`: evaluates the first derivative of the
    conformal Hubble rate $\mathrm{d}\mathcal{H} / \mathrm{d}\tau$
    * `growth_factor`: evaluates the linear matter growth factor $D_1(z)$
    * `growth_rate`: evaluates the derivative of the growth factor $f(z)$
    * `comoving_distance`: evaluates the comoving distance $\chi(z)$

    ### Main outputs

    The main objects that COFFE v3 can compute are the 2-point correlation
    function (2PCF), its multipoles, and the covariance of the multipoles.
    All of them can be computed in bulk (i.e. in a single call), using:
    * `compute_corrfunc_bulk`: for the 2PCF
    * `compute_multipoles_bulk`: for the multipoles
    * `compute_covariance_bulk`: for the covariance of the multipoles

    The output of those is an array of corresponding containers, and each
    container (`coffe.representation.Corrfunc`,
    `coffe.representation.Multipoles`, and `coffe.representation.Covariance`)
    contains the 3-dimensional coordinates and the values of the output at
    those coordinates.

    For instance, one element of the output of `compute_multipoles_bulk()` can
    be:

    ```python
    Corrfunc(
        {
            'mu': 0.99,
            'r': 95.0,
            'value': -0.004418327661669844,
            'z': 1.5
        }
    )
    ```

    Each coordinate can then be accessed using `[VARIABLE].[NAME]`, where
    `[NAME]` is one of the elements in the dictionary above (so `mu`, `r`,
    `value`, and `z`).
    The computation of the multipoles and the covariance provides analogous
    outputs.

    ### Setting coordinates for the outputs

    For the 2PCF, you can set the following coordinates:
    * `sep`: the list of comoving separations (in $\mathrm{Mpc}$) between the 2
    points in the sky
    * `mu`: the list of angles with respect to the observer
    * `z_mean`: the list of mean redshifts at which the 2PCF should be evaluated
    * `deltaz`: the list of half-widths of the redshift bins centered at `z_mean`

    For the multipoles, you can set the following coordinates:
    * `sep`: the list of comoving separations (in $\mathrm{Mpc}$) between the 2
    points in the sky
    * `l`: the list of multipoles (can be any positive integer as COFFE can
    take into account 2 populations of galaxies)
    * `z_mean`: the list of mean redshifts at which the 2PCF should be evaluated
    * `deltaz`: the list of half-widths of the redshift bins centered at `z_mean`

    For the covariance of the multipoles, along with the above for the
    multipoles, you can set the following coordinates:
    * `pixelsize`: the list of limiting sizes (resolutions) (in $\mathrm{Mpc}$)
    for each redshift bin
    * `number_density1` and `number_density2`: the list of number densities (in
    $1 / \mathrm{Mpc}^3$) for the first and second population of galaxies,
    respectively
    * `fsky`: the list of sky fractions surveyed at each redshift bin

    ### Setting contributions to outputs

    In COFFE you can control which contributions should be taken into account
    when computing the output.

    For the 2PCF and its multipoles, these are:
    * `has_density`: controls whether to take into account the density
    contribution
    * `has_rsd`: -||- the redshift-space distortion (RSD) contribution
    * `has_lensing`: -||- the magnification lensing contribution
    * `has_d1`: -||- the Doppler 1 term
    * `has_d2`: -||- the Doppler 2 term
    * `has_g1`: -||- the local gravitational term 1
    * `has_g2`: -||- the local gravitational term 2
    * `has_g3`: -||- the local gravitational term 3

    By default, COFFE computes the cross-correlation of all of the currently
    active terms + their auto correlations. This can be controlled with the
    `has_only_cross_correlations` option.

    Furthermore, one can control whether the flat-sky approximation should be
    used with the following options:
    * `has_flatsky_local`: controls whether to use the flat-sky approximation
    for local terms (such as density, RSD, and their cross-correlations)
    * `has_flatsky_local_nonlocal`: -||- for cross terms between local and
    non-local terms (such as density-lensing)
    * `has_flatsky_nonlocal`: -||- for non-local (integrated) terms (such as
    lensing-lensing)

    For the covariance of the multipoles, only the `has_density` and `has_rsd`
    contributions have an effect, as all of the other effects are not
    implemented.

    By default, COFFE computes all of the terms which contribute to the
    covariance; this can be controlled with the following options:
    * `covariance_cosmic`: controls whether to take into account the cosmic
    variance autocorrelation term
    * `covariance_poisson`: -||- the Poisson noise autocorrelation term
    * `covariance_mixed`: -||- the cross correlation term between cosmic
    variance and Poisson noise

    Note that for the covariance you can control for which populations you
    would like to compute it for with the `covariance_populations` option.

    Also note that you can specify whether or not the covariance should be
    averaged (binned) over the pixelsize or not using the
    `has_binned_covariance` option.

    ### Setting precision parameters

    You can set certain precision parameters in COFFE:
    * `integration_sampling`: how many points are sampled when computing
    2-dimensional (or higher) integrals
    * `background_sampling`: how many points are sampled between redshifts 0
    and 15 to compute the background quantities
    * `k_min`: the minimum wavenumber (in $1/\mathrm{Mpc}$) used for computing the
    Fourier-Bessel transform
    * `k_max`: the maximum wavenumber (in $1/\mathrm{Mpc}$) used for computing the
    Fourier-Bessel transform
    """
    cdef ccoffe.coffe_parameters_t _parameters
    cdef ccoffe.coffe_background_t _background
    cdef ccoffe.coffe_integral_array_t _integral
    cdef ccoffe.coffe_corrfunc_array_t _corrfunc
    cdef ccoffe.coffe_multipoles_array_t _multipoles
    cdef ccoffe.coffe_covariance_array_t _covariance_multipoles
    cdef ccoffe.coffe_covariance_array_t __dummy
    cdef int _power_spectrum_flag
    # the max size allowed for setting the various *_sampling parameters
    cdef size_t _max_size

    def __cinit__(self, **kwargs):
        """
        Constructor that initializes the structures and sets default parameters.

        Examples
        -------------
        >>> import coffe
        >>> cosmo = coffe.Coffe(omega_m=0.32, h=0.7, has_lensing=True)
        """
        # disable paralellization in Cuba, leave in only the OpenMP one
        os.environ['CUBACORES'] = '0'

        self._max_size = 2147483647

        ccoffe.coffe_parse_default_parameters(&self._parameters)

        self._background.flag = 0

        self._integral.size = 0
        self._integral.array = NULL

        self._corrfunc.size = 0
        self._corrfunc.array = NULL

        self._multipoles.size = 0
        self._multipoles.array = NULL

        self._covariance_multipoles.size = 0
        self._covariance_multipoles.array = NULL

        self._power_spectrum_flag = 0

        if kwargs:
            self.set_parameters(**kwargs)


    @staticmethod
    def from_file(filename):
        """
        Constructor that initializes the structures and sets parameters from a
        file.

        Parameters
        ----------
        filename : str
            The name of the file containing the parameters.

        Returns
        -------
        cosmo : Coffe
            The Coffe object with the parameters set from the file.

        Examples
        --------
        >>> from coffe import Coffe
        >>> cosmo = Coffe.from_file("settings.cfg")
        """
        section = "default"
        contents = f"[{section}]\n" + "\n".join(
            [_.strip(";") for _ in Path(filename).read_text().split("\n")]
        )
        config = CoffeConfigParser()
        config.read_string(contents)
        config = config[section]

        cosmo = Coffe()

        # some of the other options have different names now, or have been
        # removed; for compatibility reasons, now first the legacy options are
        # read, followed by the new ones (if they exist, they override the
        # legacy ones)
        legacy_options = [
            LegacyOption("only_cross_correlations", "int", "has_only_cross_correlations"),
            LegacyOption("covariance_z_mean", "float_array", "z_mean"),
            LegacyOption("covariance_deltaz", "float_array", "deltaz"),
            LegacyOption("covariance_fsky", "float_array", "fsky"),
            LegacyOption("covariance_pixelsize", "float_array", "pixelsize"),
            LegacyOption("covariance_step_size", "float"),
            LegacyOption("covariance_zmin", "float"),
            LegacyOption("covariance_zmax", "float"),
            LegacyOption("covariance_minimum_separation", "float"),
            LegacyOption("covariance_window", "int", "has_binned_covariance"),
            LegacyOption(
                "covariance_density", "float_array", ["number_density1", "number_density2"]
            ),
            LegacyOption("covariance_integration_bins", "int"),
            LegacyOption("multipoles", "int_array", "l"),
            LegacyOption("flatsky_local", "int", "has_flatsky_local"),
            LegacyOption("flatsky_local_nonlocal", "int", "has_flatsky_local_nonlocal"),
            LegacyOption("flatsky_nonlocal", "int", "has_flatsky_nonlocal"),
        ]

        for legacy_option in legacy_options:
            if legacy_option.name in config:
                if legacy_option.replaced_by:
                    warnings.warn(
                        f"Option '{legacy_option.name}' is no longer supported; "
                        f"it may be replaced by '{legacy_option.replaced_by}' (if found).",
                        DeprecationWarning,
                    )

                    if not isinstance(legacy_option.replaced_by, list):
                        setattr(
                            cosmo,
                            legacy_option.replaced_by,
                            getattr(config, f"get{legacy_option.kind}")(legacy_option.name),
                        )
                    else:
                        for item in legacy_option.replaced_by:
                            setattr(
                                cosmo,
                                item,
                                getattr(config, f"get{legacy_option.kind}")(legacy_option.name),
                            )
                else:
                    warnings.warn(
                        f"Option '{legacy_option.name}' is no longer supported and will be ignored.",
                        DeprecationWarning,
                    )

        if "correlation_contributions" in config:
            warnings.warn(
                f"Option 'correlation_contributions' is no longer supported; "
                f"it may be replaced by 'has_density', 'has_rsd', and 'has_lensing' (if found).",
                DeprecationWarning,
            )
            contributions = config.getstr_array("correlation_contributions")
            if "den" in contributions:
                cosmo.has_density = True
            if "rsd" in contributions:
                cosmo.has_rsd = True
            if "len" in contributions:
                cosmo.has_lensing = True

        options_float_array = [
            "sep",
            "z_mean",
            "deltaz",
            "mu",
            "number_density1",
            "number_density2",
            "pixelsize",
            "fsky",
        ]
        for option in options_float_array:
            if option in config:
                setattr(cosmo, option, config.getfloat_array(option))

        options_int_array = ["l"]
        for option in options_int_array:
            if option in config:
                setattr(cosmo, option, config.getint_array(option))

        options_float = [
            "omega_m",
            "omega_baryon",
            "omega_gamma",
            "w0",
            "wa",
            "h",
            "sigma8",
            "n_s",
            "T_cmb",
            "A_s",
            "k_per_decade_for_bao",
            "k_per_decade_for_pk",
            "start_large_k_at_tau_c_over_tau_h",
            "l_max_g",
            "l_max_ur",
            "tol_perturb_integration",
            "radiation_streaming_trigger_tau_over_tau_k",
            "ur_fluid_trigger_tau_over_tau_k",
            "N_ur",
            "m_ncdm",
            "k_min",
            "k_max",
            "YHe",
        ]
        for option in options_float:
            if option in config:
                setattr(cosmo, option, config.getfloat(option))

        options_int = [
            "has_density",
            "has_lensing",
            "has_rsd",
            "integration_sampling",
            "background_sampling",
            "has_binned_covariance",
            "has_flatsky_local",
            "has_flatsky_nonlocal",
            "has_flatsky_local_nonlocal",
            "covariance_cosmic",
            "covariance_mixed",
            "covariance_poisson",
            "N_ncdm",
        ]
        for option in options_int:
            if option in config:
                setattr(cosmo, option, config.getint(option))

        # the bias options are a bit special
        biases = ["galaxy", "magnification", "evolution"]
        populations = [1, 2]
        for bias in biases:
            for population in populations:
                if (
                    f"read_{bias}_bias{population}" in config
                    and config.getint(f"read_{bias}_bias{population}") == 1
                ):
                    getattr(cosmo, f"set_{bias}_bias{population}")(
                        *np.transpose(
                            np.loadtxt(
                                config.get(f"input_{bias}_bias{population}").strip('" ')
                            )
                        )
                    )
                elif f"input_{bias}_bias{population}" in config:
                    getattr(cosmo, f"set_{bias}_bias{population}")(
                        np.linspace(0, 10, 100),
                        [config.getfloat(f"{bias}_bias{population}")] * 100,
                    )
                elif (
                    f"{bias}_bias{population}_redshifts" in config
                    and f"{bias}_bias{population}_values" in config
                ):
                    getattr(cosmo, f"set_{bias}_bias{population}")(
                        config.getfloat_array(f"{bias}_bias{population}_redshifts"),
                        config.getfloat_array(f"{bias}_bias{population}_values"),
                    )

        # the galaxy bias can have 4 populations
        for population in [3, 4]:
            if (
                f"read_galaxy_bias{population}" in config
                and config.getint(f"read_galaxy_bias{population}") == 1
            ):
                getattr(cosmo, f"set_galaxy_bias{population}")(
                    *np.transpose(
                        np.loadtxt(
                            config.get(f"input_galaxy_bias{population}").strip('" ')
                        )
                    )
                )
            elif f"input_galaxy_bias{population}" in config:
                getattr(cosmo, f"set_galaxy_bias{population}")(
                    np.linspace(0, 10, 100),
                    [config.getfloat(f"galaxy_bias{population}")] * 100,
                )
            elif (
                f"galaxy_bias{population}_redshifts" in config
                and f"galaxy_bias{population}_values" in config
            ):
                getattr(cosmo, f"set_galaxy_bias{population}")(
                    config.getfloat_array(f"galaxy_bias{population}_redshifts"),
                    config.getfloat_array(f"galaxy_bias{population}_values"),
                )

        return cosmo


    def to_file(self, filename):
        """Write the cosmology to a file.

        Parameters
        ----------
        filename : str
            The name of the file.
        """
        with open(filename, "w", encoding="utf-8") as f:
            for key, value in self.parameters.items():
                if isinstance(value, bool):
                    value = int(value)
                if isinstance(value, (list, tuple, np.ndarray)):
                    value = list(value)
                f.write(f"{key} = {value}\n")

            # the various biases are handled separately

            # -------------------------------------------------------------------------------
            # galaxy bias
            # -------------------------------------------------------------------------------
            f.write(f"galaxy_bias1_redshifts = {self.galaxy_bias1_redshifts().tolist()}\n")
            f.write(f"galaxy_bias1_values = {self.galaxy_bias1_values().tolist()}\n")

            f.write(f"galaxy_bias2_redshifts = {self.galaxy_bias2_redshifts().tolist()}\n")
            f.write(f"galaxy_bias2_values = {self.galaxy_bias2_values().tolist()}\n")

            f.write(f"galaxy_bias3_redshifts = {self.galaxy_bias3_redshifts().tolist()}\n")
            f.write(f"galaxy_bias3_values = {self.galaxy_bias3_values().tolist()}\n")

            f.write(f"galaxy_bias4_redshifts = {self.galaxy_bias4_redshifts().tolist()}\n")
            f.write(f"galaxy_bias4_values = {self.galaxy_bias4_values().tolist()}\n")

            # -------------------------------------------------------------------------------
            # magnification bias
            # -------------------------------------------------------------------------------
            f.write(f"magnification_bias1_redshifts = {self.magnification_bias1_redshifts().tolist()}\n")
            f.write(f"magnification_bias1_values = {self.magnification_bias1_values().tolist()}\n")

            f.write(f"magnification_bias2_redshifts = {self.magnification_bias2_redshifts().tolist()}\n")
            f.write(f"magnification_bias2_values = {self.magnification_bias2_values().tolist()}\n")

            # -------------------------------------------------------------------------------
            # evolution bias
            # -------------------------------------------------------------------------------
            f.write(f"evolution_bias1_redshifts = {self.evolution_bias1_redshifts().tolist()}\n")
            f.write(f"evolution_bias1_values = {self.evolution_bias1_values().tolist()}\n")

            f.write(f"evolution_bias2_redshifts = {self.evolution_bias2_redshifts().tolist()}\n")
            f.write(f"evolution_bias2_values = {self.evolution_bias2_values().tolist()}\n")


    def galaxy_bias1_redshifts(self):
        """Return the redshifts of the galaxy bias 1 spline."""
        size = get_spline_size(&self._parameters.galaxy_bias1)
        xmin = get_spline_min(&self._parameters.galaxy_bias1)
        xmax = get_spline_max(&self._parameters.galaxy_bias1)

        return np.linspace(xmin, xmax, size)

    def galaxy_bias1_values(self):
        """Return the values of the galaxy bias 1 spline."""
        return np.array([self.galaxy_bias1(_) for _ in self.galaxy_bias1_redshifts()])


    def galaxy_bias2_redshifts(self):
        """Return the redshifts of the galaxy bias 2 spline."""
        size = get_spline_size(&self._parameters.galaxy_bias2)
        xmin = get_spline_min(&self._parameters.galaxy_bias2)
        xmax = get_spline_max(&self._parameters.galaxy_bias2)

        return np.linspace(xmin, xmax, size)

    def galaxy_bias2_values(self):
        """Return the values of the galaxy bias 2 spline."""
        return np.array([self.galaxy_bias2(_) for _ in self.galaxy_bias2_redshifts()])

    def galaxy_bias3_redshifts(self):
        """Return the redshifts of the galaxy bias 3 spline."""
        size = get_spline_size(&self._parameters.galaxy_bias3)
        xmin = get_spline_min(&self._parameters.galaxy_bias3)
        xmax = get_spline_max(&self._parameters.galaxy_bias3)

        return np.linspace(xmin, xmax, size)

    def galaxy_bias3_values(self):
        """Return the values of the galaxy bias 3 spline."""
        return np.array([self.galaxy_bias3(_) for _ in self.galaxy_bias3_redshifts()])

    def galaxy_bias4_redshifts(self):
        """Return the redshifts of the galaxy bias 4 spline."""
        size = get_spline_size(&self._parameters.galaxy_bias4)
        xmin = get_spline_min(&self._parameters.galaxy_bias4)
        xmax = get_spline_max(&self._parameters.galaxy_bias4)

        return np.linspace(xmin, xmax, size)

    def galaxy_bias4_values(self):
        """Return the values of the galaxy bias 4 spline."""
        return np.array([self.galaxy_bias4(_) for _ in self.galaxy_bias4_redshifts()])


    def magnification_bias1_redshifts(self):
        """Return the redshifts of the magnification bias 1 spline."""
        size = get_spline_size(&self._parameters.magnification_bias1)
        xmin = get_spline_min(&self._parameters.magnification_bias1)
        xmax = get_spline_max(&self._parameters.magnification_bias1)

        return np.linspace(xmin, xmax, size)

    def magnification_bias1_values(self):
        """Return the values of the magnification bias 1 spline."""
        return np.array([self.magnification_bias1(_) for _ in self.magnification_bias1_redshifts()])


    def magnification_bias2_redshifts(self):
        """Return the redshifts of the magnification bias 2 spline."""
        size = get_spline_size(&self._parameters.magnification_bias2)
        xmin = get_spline_min(&self._parameters.magnification_bias2)
        xmax = get_spline_max(&self._parameters.magnification_bias2)

        return np.linspace(xmin, xmax, size)

    def magnification_bias2_values(self):
        """Return the values of the magnification bias 2 spline."""
        return np.array([self.magnification_bias2(_) for _ in self.magnification_bias2_redshifts()])


    def evolution_bias1_redshifts(self):
        """Return the redshifts of the evolution bias 1 spline."""
        size = get_spline_size(&self._parameters.evolution_bias1)
        xmin = get_spline_min(&self._parameters.evolution_bias1)
        xmax = get_spline_max(&self._parameters.evolution_bias1)

        return np.linspace(xmin, xmax, size)

    def evolution_bias1_values(self):
        """Return the values of the evolution bias 1 spline."""
        return np.array([self.evolution_bias1(_) for _ in self.evolution_bias1_redshifts()])


    def evolution_bias2_redshifts(self):
        """Return the redshifts of the evolution bias 2 spline."""
        size = get_spline_size(&self._parameters.evolution_bias2)
        xmin = get_spline_min(&self._parameters.evolution_bias2)
        xmax = get_spline_max(&self._parameters.evolution_bias2)

        return np.linspace(xmin, xmax, size)

    def evolution_bias2_values(self):
        """Return the values of the evolution bias 2 spline."""
        return np.array([self.evolution_bias2(_) for _ in self.evolution_bias2_redshifts()])


    def _free_background(self):
        ccoffe.coffe_background_free(&self._background)


    def _free_power_spectrum(self):
        ccoffe.coffe_free_spline(&self._parameters.power_spectrum)
        ccoffe.coffe_free_spline(&self._parameters.power_spectrum_norm)
        ccoffe.coffe_free_spline2d(&self._parameters.power_spectrum2d)
        ccoffe.coffe_free_spline2d(&self._parameters.power_spectrum2d_norm)
        self._power_spectrum_flag = 0


    def _free_integrals(self):
        ccoffe.coffe_integrals_free(&self._integral)


    def _free_corrfunc(self):
        ccoffe.coffe_corrfunc_free(&self._corrfunc)


    def _free_multipoles(self):
        ccoffe.coffe_multipoles_free(&self._multipoles)


    def _free_covariance_multipoles(self):
        ccoffe.coffe_covariance_free(&self._covariance_multipoles)


    def _free_except_parameters(self):
        self._free_background()
        self._free_power_spectrum()
        self._free_integrals()
        self._free_corrfunc()
        self._free_multipoles()
        self._free_covariance_multipoles()


    def __dealloc__(self):
        self._free_except_parameters()
        ccoffe.coffe_parameters_free(&self._parameters)


    @property
    def verbose(self):
        return bool(self._parameters.verbose)

    @verbose.setter
    def verbose(self, value):
        self._parameters.verbose = int(bool(value))


    @property
    def use_little_omega(self):
        try:
            return bool(int(os.environ.get('COFFE_USE_LITTLE_OMEGA')))
        except:
            return False


    @property
    def _coeff(self):
        return self.h**2 if self.use_little_omega else 1


    @property
    def has_cuba(self):
        return bool(self._parameters.has_cuba)

    @property
    def has_class(self):
        return bool(self._parameters.has_class)

    @property
    def big_omega_m(self):
        return self._parameters.Omega0_m

    @property
    def big_omega_de(self):
        return self._parameters.Omega0_de

    @property
    def big_omega_baryon(self):
        return self._parameters.Omega0_baryon

    @property
    def big_omega_cdm(self):
        return self._parameters.Omega0_cdm

    @property
    def big_omega_nu(self):
        return self._parameters.Omega0_nu

    @property
    def big_omega_gamma(self):
        return self._parameters.Omega0_gamma


    @staticmethod
    def _check_omegas(*omegas):
        if np.sum(np.array([_ for _ in omegas])) > 1:
            raise ValueError(
                f'Cannot satisfy the condition Omega_de >= 0'
            )


    def _check_coords_corrfunc(self):
        # check there is actually something to compute
        if not (len(self.mu) and len(self.sep) and len(self.z_mean)):
            raise ValueError(
                f'At least one of the following parameters is empty: mu ({self.mu}), sep ({self.sep}), z_mean ({self.z_mean})'
            )


    def _check_coords_multipoles(self):
        # check there is actually something to compute
        if not (len(self.l) and len(self.sep) and len(self.z_mean)):
            raise ValueError(
                f'At least one of the following parameters is empty: l ({self.l}), sep ({self.sep}), z_mean ({self.z_mean})'
            )


    def _check_contributions(self):
        if not any(
            [
                self.has_density,
                self.has_rsd,
                self.has_lensing,
                self.has_d1,
                self.has_d2,
                self.has_g1,
                self.has_g2,
                self.has_g3,
                self.has_g4,
                self.has_g5,
            ]
        ):
            raise ValueError(
                "No contributions specified, you need to specify at least one of: "
                "'has_density', 'has_rsd', 'has_lensing, "
                "'has_d1', 'has_d2', 'has_g1', 'has_g2', 'has_g3', "
                "'has_g4', 'has_g5'"
            )
        if any([self.has_d2, self.has_g1, self.has_g2, self.has_g3, self.has_g4, self.has_g5]):
            self._parameters.divergent = 1


    def _balance_content(self):
        """
        Internal function that balances the energy content of dark energy so it all adds up to 1.
        Called automatically when any setter for the various omegas is called.
        """
        self._parameters.Omega0_de = 1 - self._parameters.Omega0_m - self._parameters.Omega0_gamma


    @property
    def has_only_cross_correlations(self):
        """
        Whether or not we consider only cross-correlations.
        """
        return bool(self._parameters.only_cross_correlations)

    @has_only_cross_correlations.setter
    def has_only_cross_correlations(self, value):
        self._parameters.only_cross_correlations = int(bool(value))
        self._free_corrfunc()
        self._free_multipoles()


    @property
    def parameters(self):
        """
        Returns the current writable parameters as a dictionary (can be
        re-used for `set_parameters`).
        """
        properties = {
            key : getattr(self, key) \
            for key in dir(self.__class__) \
            if key != 'parameters' and hasattr(getattr(self.__class__, key), '__set__')
        }
        writable = {}
        for key in properties:
            try:
                setattr(self, key, getattr(self, key))
                writable[key] = getattr(self, key)
            except AttributeError:
                pass
        return writable


    def set_parameters(self, **value):
        """
        Bulk setter of parameters.
        """
        for key in value:
            if not hasattr(self, key):
                raise AttributeError(
                    f'The parameter {key} is not a valid parameter in COFFE.'
                )

        for key in value:
            setattr(self, key, value[key])


    @property
    def omega_cdm(self):
        """
        Fraction of cold dark matter today.
        """
        return self._parameters.Omega0_cdm * self._coeff

    @omega_cdm.setter
    def omega_cdm(self, value):
        _check_parameter('omega_cdm', value, (int, float), -100, 100)
        if not np.allclose(value, self.omega_cdm):
            # we set the value, rebalance the Omega budget, and free memory
            self._check_omegas(value, self.omega_baryon, self.omega_gamma, self.omega_nu)
            self._parameters.Omega0_cdm = value / self._coeff
            self._parameters.Omega0_m = self._parameters.Omega0_cdm + self._parameters.Omega0_baryon
            self._balance_content()
            self._free_except_parameters()


    @property
    def omega_nu(self):
        """
        Returns the energy density fraction of massive neutrinos.
        """
        return self._parameters.Omega0_nu * self._coeff


    @property
    def omega_m(self):
        """
        Fraction of total matter (CDM + baryons) today.
        """
        return self._parameters.Omega0_m * self._coeff

    @omega_m.setter
    def omega_m(self, value):
        _check_parameter('omega_m', value, (int, float), -100, 100)
        if not np.allclose(value, self.omega_m):
            # we set the value, rebalance the Omega budget, and free memory
            self._check_omegas(value, self.omega_gamma)
            self._parameters.Omega0_m = value / self._coeff
            # for consistency with TotallySAF
            self._parameters.Omega0_cdm = self._parameters.Omega0_m - self._parameters.Omega0_baryon - self._parameters.Omega0_nu
            self._balance_content()
            self._free_except_parameters()


    @property
    def omega_baryon(self):
        """
        Fraction of baryonic matter today.
        """
        return self._parameters.Omega0_baryon * self._coeff

    @omega_baryon.setter
    def omega_baryon(self, value):
        _check_parameter('omega_baryon', value, (int, float), -100, 100)
        if not np.allclose(value, self.omega_baryon):
            # we set the value, rebalance the Omega budget, and free memory
            self._check_omegas(value, self.omega_cdm, self.omega_gamma)
            self._parameters.Omega0_baryon = value / self._coeff
            self._parameters.Omega0_cdm = self._parameters.Omega0_m - self._parameters.Omega0_baryon - self._parameters.Omega0_nu
            self._balance_content()
            self._free_except_parameters()


    @property
    def omega_gamma(self):
        """
        Fraction of relativistic species today.
        """
        return self._parameters.Omega0_gamma * self._coeff

    @omega_gamma.setter
    def omega_gamma(self, value):
        _check_parameter('omega_gamma', value, (int, float), -100, 100)
        if not np.allclose(value, self.omega_gamma):
            self._check_omegas(value, self.omega_m)
            self._parameters.Omega0_gamma = value / self._coeff
            self._balance_content()
            self._free_except_parameters()


    @property
    def omega_de(self):
        """
        Dark energy fraction today.
        """
        return self._parameters.Omega0_de * self._coeff


    @property
    def h(self):
        """
        The 'little h' parameter (reduced Hubble rate).
        """
        return self._parameters.h

    @h.setter
    def h(self, value):
        _check_parameter('h', value, (int, float), 0, 100)
        if not np.allclose(value, self.h):
            # setting the old value equal to the new value of h
            self._parameters.h = value
            self._balance_content()
            self._free_except_parameters()


    @property
    def w0(self):
        """
        Dark energy is parametrized by the equation of state:

        .. math::
            w(z) = w_0 + (1 - a) w_a

        """
        return self._parameters.w0

    @w0.setter
    def w0(self, value):
        _check_parameter('w0', value, (int, float), -100, 100)
        if not np.allclose(value, self.w0):
            self._parameters.w0 = value
            self._free_except_parameters()


    @property
    def k_per_decade_for_bao(self):
        return self._parameters.class_precision.k_per_decade_for_bao

    @k_per_decade_for_bao.setter
    def k_per_decade_for_bao(self, value):
        self._parameters.class_precision.k_per_decade_for_bao = value
        self._free_except_parameters()


    @property
    def k_per_decade_for_pk(self):
        return self._parameters.class_precision.k_per_decade_for_pk

    @k_per_decade_for_pk.setter
    def k_per_decade_for_pk(self, value):
        self._parameters.class_precision.k_per_decade_for_pk = value
        self._free_except_parameters()


    @property
    def start_large_k_at_tau_h_over_tau_k(self):
        return self._parameters.class_precision.start_large_k_at_tau_h_over_tau_k

    @start_large_k_at_tau_h_over_tau_k.setter
    def start_large_k_at_tau_h_over_tau_k(self, value):
        self._parameters.class_precision.start_large_k_at_tau_h_over_tau_k = value
        self._free_except_parameters()


    @property
    def l_max_g(self):
        return self._parameters.class_precision.l_max_g

    @l_max_g.setter
    def l_max_g(self, value):
        self._parameters.class_precision.l_max_g = value
        self._free_except_parameters()


    @property
    def l_max_ur(self):
        return self._parameters.class_precision.l_max_ur

    @l_max_ur.setter
    def l_max_ur(self, value):
        self._parameters.class_precision.l_max_ur = value
        self._free_except_parameters()


    @property
    def tol_perturb_integration(self):
        return self._parameters.class_precision.tol_perturb_integration

    @tol_perturb_integration.setter
    def tol_perturb_integration(self, value):
        self._parameters.class_precision.tol_perturb_integration = value
        self._free_except_parameters()


    @property
    def radiation_streaming_trigger_tau_over_tau_k(self):
        return self._parameters.class_precision.radiation_streaming_trigger_tau_over_tau_k

    @radiation_streaming_trigger_tau_over_tau_k.setter
    def radiation_streaming_trigger_tau_over_tau_k(self, value):
        self._parameters.class_precision.radiation_streaming_trigger_tau_over_tau_k = value
        self._free_except_parameters()


    @property
    def ur_fluid_trigger_tau_over_tau_k(self):
        return self._parameters.class_precision.ur_fluid_trigger_tau_over_tau_k

    @ur_fluid_trigger_tau_over_tau_k.setter
    def ur_fluid_trigger_tau_over_tau_k(self, value):
        self._parameters.class_precision.ur_fluid_trigger_tau_over_tau_k = value
        self._free_except_parameters()


    @property
    def wa(self):
        """
        Dark energy is parametrized by the equation of state:

        .. math::
            w(z) = w_0 + (1 - a) w_a

        """
        return self._parameters.wa

    @wa.setter
    def wa(self, value):
        _check_parameter('wa', value, (int, float), -100, 100)
        if not np.allclose(value, self.wa):
            self._parameters.wa = value
            self._free_except_parameters()


    @property
    def n_s(self):
        """
        Spectral index (tilt) of the primordial power spectrum.
        """
        return self._parameters.n_s

    @n_s.setter
    def n_s(self, value):
        _check_parameter('n_s', value, (int, float), -100, 100)
        if not np.allclose(value, self.n_s):
            self._parameters.n_s = value
            self._free_except_parameters()


    @property
    def sigma8(self):
        r"""
        Windowed density fluctuation at $r = 8\ \mathrm{Mpc}/h$
        """
        return self._parameters.sigma8

    @sigma8.setter
    def sigma8(self, value):
        _check_parameter('sigma8', value, (int, float), -100, 100)
        if not np.allclose(value, self.sigma8):
            self._parameters.sigma8 = value
            self._free_power_spectrum()
            self._free_integrals()
            self._free_corrfunc()
            self._free_multipoles()
            self._free_covariance_multipoles()


    @property
    def A_s(self):
        """
        The amplitude of the primordial power spectrum
        """
        return self._parameters.A_s

    @A_s.setter
    def A_s(self, value):
        _check_parameter('A_s', value, (int, float), -100, 100)
        self._parameters.A_s = value
        self._free_power_spectrum()
        self._free_integrals()
        self._free_corrfunc()
        self._free_multipoles()
        self._free_covariance_multipoles()


    @property
    def T_cmb(self):
        """
        The average temperature of the CMB.
        """
        return self._parameters.T_cmb

    @T_cmb.setter
    def T_cmb(self, value):
        _check_parameter('T_cmb', value, (int, float), -100, 100)
        if not np.allclose(value, self.T_cmb):
            self._parameters.T_cmb = value
            self._free_except_parameters()


    @property
    def N_ur(self):
        """
        The number of ultra-relativistic species.
        """
        return self._parameters.N_ur

    @N_ur.setter
    def N_ur(self, value):
        _check_parameter('N_ur', value, (int, float), -100, 100)
        if not np.allclose(value, self.N_ur):
            self._parameters.N_ur = value
            self._free_except_parameters()


    @property
    def m_ncdm(self):
        """
        The sum of masses of non-CDM species (mostly for neutrinos), in units
        of eV.
        """
        return self._parameters.m_ncdm

    @m_ncdm.setter
    def m_ncdm(self, value):
        _check_parameter('m_ncdm', value, (int, float), -100, 100)
        if not np.allclose(value, self.m_ncdm):
            self._parameters.m_ncdm = value
            # see eq. (19) of https://arxiv.org/abs/1212.6154
            # note that little omega(nu) doesn't change if we change h, but big
            # omega(nu) does!
            self._parameters.Omega0_nu = value / 93.14 / self.h / self.h
            self._parameters.Omega0_cdm = self._parameters.Omega0_m - self._parameters.Omega0_baryon - self._parameters.Omega0_nu
            self._balance_content()
            self._free_except_parameters()


    @property
    def YHe(self):
        """
        The primordial helium fraction
        """
        return self._parameters.YHe

    @YHe.setter
    def YHe(self, value):
        _check_parameter('YHe', value, (int, float), -100, 100)
        if not np.allclose(value, self.YHe):
            self._parameters.YHe = value
            self._free_except_parameters()


    @property
    def N_ncdm(self):
        """
        The number of (massive!) non-CDM species.
        """
        return self._parameters.N_ncdm

    @N_ncdm.setter
    def N_ncdm(self, value):
        _check_parameter('N_ncdm', value, int, -100, 100)
        if not np.allclose(value, self.N_ncdm):
            self._parameters.N_ncdm = value
            self._free_except_parameters()


    @property
    def sep(self):
        r"""
        Returns the list of separations (in $\mathrm{Mpc}$) for which the
        2PCF/multipoles/covariance of multipoles should be computed.
        """
        return np.array(
            [self._parameters.sep[i] for i in range(self._parameters.sep_len)]
        )


    @sep.setter
    def sep(self, value):
        try:
            _ = iter(value)
        except TypeError as err:
            raise TypeError(f'The value {value} is not iterable') from err

        try:
            temp = [float(_) for _ in value]
        except TypeError as err:
            raise TypeError('Cannot convert all values to floats') from err

        if not all(_> 0 and _ < 10000 for _ in temp):
            raise ValueError

        if self._parameters.sep_len:
            free(self._parameters.sep)

        self._parameters.sep_len = len(temp)
        self._parameters.sep = <double *> malloc(sizeof(double) * len(temp))
        for i in range(self._parameters.sep_len):
            self._parameters.sep[i] = temp[i]
        self._free_corrfunc()
        self._free_multipoles()


    @property
    def mu(self):
        r"""
        The list of angles mu for which the 2PCF should be computed.
        it is defined as:

        .. math::
            \mu = \frac{r_\parallel}{r} = \frac{\chi_2 - \chi_1}{r}

        """
        return np.array(
            [self._parameters.mu[i] for i in range(self._parameters.mu_len)]
        )

    @mu.setter
    def mu(self, value):
        try:
            _ = iter(value)
        except TypeError as err:
            raise TypeError(f'The value {value} is not iterable') from err

        try:
            temp = [float(_) for _ in value]
        except TypeError as err:
            raise TypeError('Cannot convert all values to floats') from err

        if not all(_> -1 and _ < 1 for _ in temp):
            raise ValueError

        if self._parameters.mu_len:
            free(self._parameters.mu)

        self._parameters.mu_len = len(temp)
        self._parameters.mu = <double *> malloc(sizeof(double) * len(temp))
        for i in range(self._parameters.mu_len):
            self._parameters.mu[i] = temp[i]
        # we need to re-compute it since we changed the values
        # TODO make it so that if we give identical values, we don't recompute it
        self._free_corrfunc()


    def integral(self, *, r : float, n : int, l : int):
        r"""
        Returns one of the Fourier Bessel integrals (only integer arguments for now).
        In essence, computes the following integral:

        .. math::
            \frac{1}{2 \pi^2} \int\limits_0^\infty dk\, k^2\, P(k, z = 0) \frac{j_\ell(k r)}{(k r)^n}

        Note that if $n > \ell$, it returns the above multiplied by $r^{n -
        \ell}$ instead (see
        [arXiv:1806.11090](https://arxiv.org/abs/1806.11090), section 4.2 for
        further explanation).

        Parameters
        ----------
        r : float
            the separation (in $\mathrm{Mpc}$) for which one wants to compute
            the integral

        n : int
            the power of the denominator

        l : int
            the order of the Bessel function

        Returns
        -------
        float
            the value of the integral
        """
        _check_parameter('r', r, (int, float), 0, 25000)
        _check_parameter('n', n, int, -2, 4)
        _check_parameter('l', l, int, 0, 4)

        if not self._background.flag:
            self._background_init()

        if not self._power_spectrum_flag:
            self._power_spectrum_init()

        if not self._integral.size:
            self._integrals_init()

        cdef ccoffe.coffe_integral_t *result = ccoffe.coffe_find_integral(
            &self._integral,
            n,
            l,
            ccoffe.COFFE_INTEGER,
            ccoffe.COFFE_INTEGER
        )

        if result == NULL:
            raise ValueError(
                f'Cannot find integral with n = {n}, l = {l}'
            )

        return ccoffe.coffe_interp_spline(&result.result, r)


    def galaxy_bias1(self, z : float):
        """
        Evaluates the galaxy bias of the first population of tracers at some
        redshift.

        Returns
        -------
        float
        """
        return evaluate_spline(&self._parameters.galaxy_bias1, z)


    def set_galaxy_bias1(self, x_sampling : List[float], y_sampling : List[float]):
        """
        Sets the value of the galaxy bias for the first population of tracers.

        Parameters
        ----------
        x_sampling : array_like of float
            the redshifts at which the bias is sampled

        y_sampling : array_like of float
            the values of the bias at the redshifts given in x_sampling
        """
        set_spline(
            &self._parameters.galaxy_bias1,
            x_sampling, y_sampling,
            self._parameters.interp_method
        )
        self._free_corrfunc()
        self._free_multipoles()


    def galaxy_bias2(self, z : float):
        """
        Evaluates the galaxy bias of the second population of tracers at some
        redshift.

        Returns
        -------
        float
        """
        return evaluate_spline(&self._parameters.galaxy_bias2, z)


    def set_galaxy_bias2(self, x_sampling : List[float], y_sampling : List[float]):
        """
        Sets the value of the galaxy bias for the second population of tracers.

        Parameters
        ----------
        x_sampling : array_like of float
            the redshifts at which the bias is sampled

        y_sampling : array_like of float
            the values of the bias at the redshifts given in x_sampling
        """
        set_spline(
            &self._parameters.galaxy_bias2,
            x_sampling, y_sampling,
            self._parameters.interp_method
        )
        self._free_corrfunc()
        self._free_multipoles()


    def galaxy_bias3(self, z : float):
        """
        Evaluates the galaxy bias of the first population of tracers at some
        redshift.
        """
        return evaluate_spline(&self._parameters.galaxy_bias3, z)


    def set_galaxy_bias3(self, x_sampling : List[float], y_sampling : List[float]):
        """
        Sets the value of the galaxy bias for the third population of tracers
        (only relevant for covariance).
        """
        set_spline(
            &self._parameters.galaxy_bias3,
            x_sampling, y_sampling,
            self._parameters.interp_method
        )
        self._free_corrfunc()
        self._free_multipoles()


    def galaxy_bias4(self, z : float):
        """
        Evaluates the galaxy bias of the first population of tracers at some
        redshift.
        """
        return evaluate_spline(&self._parameters.galaxy_bias4, z)


    def set_galaxy_bias4(self, x_sampling : List[float], y_sampling : List[float]):
        """
        Sets the value of the galaxy bias for the fourth population of tracers
        (only relevant for covariance).
        """
        set_spline(
            &self._parameters.galaxy_bias4,
            x_sampling, y_sampling,
            self._parameters.interp_method
        )
        self._free_corrfunc()
        self._free_multipoles()


    def magnification_bias1(self, z : float):
        """
        Evaluates the magnification bias of the first population at some redshift.

        Returns
        -------
        float
        """
        return evaluate_spline(&self._parameters.magnification_bias1, z)


    def set_magnification_bias1(self, x_sampling : List[float], y_sampling : List[float]):
        """
        Sets the value of the magnification bias for the first population of
        tracers.

        Parameters
        ----------
        x_sampling : array_like of float
            the redshifts at which the bias is sampled

        y_sampling : array_like of float
            the values of the bias at the redshifts given in x_sampling
        """
        set_spline(
            &self._parameters.magnification_bias1,
            x_sampling, y_sampling,
            self._parameters.interp_method
        )
        self._free_corrfunc()
        self._free_multipoles()


    def magnification_bias2(self, z : float):
        """
        Evaluates the magnification bias of the second population at some redshift.

        Returns
        -------
        float
        """
        return evaluate_spline(&self._parameters.magnification_bias2, z)


    def set_magnification_bias2(self, x_sampling : List[float], y_sampling : List[float]):
        """
        Sets the value of the magnification bias of the second population of
        tracers.

        Parameters
        ----------
        x_sampling : array_like of float
            the redshifts at which the bias is sampled

        y_sampling : array_like of float
            the values of the bias at the redshifts given in x_sampling
        """
        set_spline(
            &self._parameters.magnification_bias2,
            x_sampling, y_sampling,
            self._parameters.interp_method
        )
        self._free_corrfunc()
        self._free_multipoles()


    def evolution_bias1(self, z : float):
        """
        Evaluates the evolution bias of the first population at some redshift.
        """
        return evaluate_spline(&self._parameters.evolution_bias1, z)


    def set_evolution_bias1(self, x_sampling : List[float], y_sampling : List[float]):
        """
        Sets the value of the evolution bias for the first population of
        tracers.
        """
        set_spline(
            &self._parameters.evolution_bias1,
            x_sampling, y_sampling,
            self._parameters.interp_method
        )
        self._free_corrfunc()
        self._free_multipoles()


    def evolution_bias2(self, z : float):
        """
        Evaluates the evolution bias of the second population at some redshift.
        """
        return evaluate_spline(&self._parameters.evolution_bias2, z)


    def set_evolution_bias2(self, x_sampling : List[float], y_sampling : List[float]):
        """
        Sets the value of the evolution bias of the second population of
        tracers.
        """
        set_spline(
            &self._parameters.evolution_bias2,
            x_sampling, y_sampling,
            self._parameters.interp_method
        )
        self._free_corrfunc()
        self._free_multipoles()


    @property
    def l(self):
        """
        Returns the multipole moments for which the multipoles should be computed.
        """
        return np.array(
            [self._parameters.multipole_values[i] for i in range(self._parameters.multipole_values_len)]
        )


    @l.setter
    def l(self, value):
        try:
            _ = iter(value)
        except TypeError as err:
            raise TypeError(f'The value {value} is not iterable') from err

        try:
            temp = [_.__index__() for _ in value]
        except AttributeError as err:
            raise TypeError('Cannot convert all values to int') from err

        if not all(_ >= 0 and _ < 10 for _ in temp):
            raise ValueError

        if self._parameters.multipole_values_len:
            free(self._parameters.multipole_values)

        self._parameters.multipole_values_len = len(temp)
        self._parameters.multipole_values = <int *> malloc(sizeof(int) * len(temp))
        for i in range(self._parameters.multipole_values_len):
            self._parameters.multipole_values[i] = temp[i]
        # we need to free the integrals since the flat-sky ones may depend on them
        self._free_multipoles()
        self._free_covariance_multipoles()


    @property
    def z_mean(self):
        """
        The mean redshift for which the signal (2PCF/multipoles) or covariance
        should be computed.
        """
        return np.array(
            [self._parameters.z_mean[i] for i in range(self._parameters.z_mean_len)]
        )

    @z_mean.setter
    def z_mean(self, value):
        try:
            _ = iter(value)
        except TypeError as err:
            raise TypeError(f'The value {value} is not iterable') from err

        try:
            temp = [float(_) for _ in value]
        except TypeError as err:
            raise TypeError('Cannot convert all values to floats') from err

        if not all(_>= 0 and _ < 10 for _ in temp):
            raise ValueError

        if self._parameters.z_mean_len:
            free(self._parameters.z_mean)

        self._parameters.z_mean_len = len(temp)
        self._parameters.z_mean = <double *> malloc(sizeof(double) * len(temp))
        for i in range(self._parameters.z_mean_len):
            self._parameters.z_mean[i] = temp[i]
        self._free_corrfunc()
        self._free_multipoles()
        self._free_covariance_multipoles()


    @property
    def number_density1(self):
        r"""
        Number density of first tracers (in $1/\mathrm{Mpc}^3$) at `z_mean`.
        """
        return np.array(
            [self._parameters.density1[i] for i in range(self._parameters.density1_len)]
        )

    @number_density1.setter
    def number_density1(self, value):
        try:
            _ = iter(value)
        except TypeError as err:
            raise TypeError(f'The value {value} is not iterable') from err

        try:
            temp = [float(_) for _ in value]
        except TypeError as err:
            raise TypeError('Cannot convert all values to floats') from err

        if not all(_>= 0 and _ < 10 for _ in temp):
            raise ValueError

        if self._parameters.density1_len:
            free(self._parameters.density1)

        self._parameters.density1_len = len(temp)
        self._parameters.density1 = <double *> malloc(sizeof(double) * len(temp))
        for i in range(self._parameters.density1_len):
            self._parameters.density1[i] = temp[i]
        self._free_corrfunc()
        self._free_multipoles()
        self._free_covariance_multipoles()


    @property
    def number_density2(self):
        """
        Number density of second tracers (in $1/\mathrm{Mpc}^3$) at `z_mean`.
        """
        return np.array(
            [self._parameters.density2[i] for i in range(self._parameters.density2_len)]
        )

    @number_density2.setter
    def number_density2(self, value):
        try:
            _ = iter(value)
        except TypeError as err:
            raise TypeError(f'The value {value} is not iterable') from err

        try:
            temp = [float(_) for _ in value]
        except TypeError as err:
            raise TypeError('Cannot convert all values to floats') from err

        if not all(_>= 0 and _ < 10 for _ in temp):
            raise ValueError

        if self._parameters.density2_len:
            free(self._parameters.density2)

        self._parameters.density2_len = len(temp)
        self._parameters.density2 = <double *> malloc(sizeof(double) * len(temp))
        for i in range(self._parameters.density2_len):
            self._parameters.density2[i] = temp[i]
        self._free_corrfunc()
        self._free_multipoles()
        self._free_covariance_multipoles()


    @property
    def covariance_populations(self):
        """
        Returns the populations used for the covariance matrix.
        """
        return np.array([
            self._parameters.covariance_pop1,
            self._parameters.covariance_pop2,
            self._parameters.covariance_pop3,
            self._parameters.covariance_pop4,
        ])

    @covariance_populations.setter
    def covariance_populations(self, value):
        try:
            _ = iter(value)
        except TypeError as err:
            raise TypeError(f'The value {value} is not iterable') from err

        try:
            temp = [int(_) for _ in value]
        except TypeError as err:
            raise TypeError('Cannot convert all values to integers') from err

        if len(value) != 4:
            raise ValueError('The array for setting `covariance_populations` must have length 4')

        if any(_ not in [1, 2, 3, 4] for _ in value):
            raise ValueError('The values for `covariance_populations` must be one of: [1, 2, 3, 4]')

        self._parameters.covariance_pop1, self._parameters.covariance_pop2, self._parameters.covariance_pop3, self._parameters.covariance_pop4 = value
        self._free_covariance_multipoles()


    @property
    def covariance_cosmic(self):
        """
        Whether the CV-CV term should be taken into account when computing the
        covariance of multipoles
        """
        return bool(self._parameters.covariance_cosmic)

    @covariance_cosmic.setter
    def covariance_cosmic(self, value):
        self._parameters.covariance_cosmic = int(bool(value))
        self._free_covariance_multipoles()


    @property
    def covariance_mixed(self):
        """
        Whether the mixed (CV-Poisson + Poisson-CV) term should be taken into account when computing the
        covariance of multipoles
        """
        return bool(self._parameters.covariance_mixed)

    @covariance_mixed.setter
    def covariance_mixed(self, value):
        self._parameters.covariance_mixed = int(bool(value))
        self._free_covariance_multipoles()

    @property
    def covariance_poisson(self):
        """
        Whether the Poisson-Poisson term should be taken into account when computing the
        covariance of multipoles
        """
        return bool(self._parameters.covariance_poisson)

    @covariance_poisson.setter
    def covariance_poisson(self, value):
        self._parameters.covariance_poisson = int(bool(value))
        self._free_covariance_multipoles()


    @property
    def pixelsize(self):
        r"""
        The pixel size of the covariance (roughly the resolution of the survey)
        in $\mathrm{Mpc}$.
        """
        return np.array(
            [self._parameters.pixelsize[i] for i in range(self._parameters.pixelsize_len)]
        )

    @pixelsize.setter
    def pixelsize(self, value):
        try:
            _ = iter(value)
        except TypeError as err:
            raise TypeError(f'The value {value} is not iterable') from err

        try:
            temp = [float(_) for _ in value]
        except TypeError as err:
            raise TypeError('Cannot convert all values to floats') from err

        if not all(_> 0 and _ < 1000 for _ in temp):
            raise ValueError

        if self._parameters.pixelsize_len:
            free(self._parameters.pixelsize)

        self._parameters.pixelsize_len = len(temp)
        self._parameters.pixelsize = <double *> malloc(sizeof(double) * len(temp))
        for i in range(self._parameters.pixelsize_len):
            self._parameters.pixelsize[i] = temp[i]
        self._free_corrfunc()
        self._free_multipoles()
        self._free_covariance_multipoles()


    @property
    def fsky(self):
        """
        The sky fraction covered by the survey.
        Must be given as a list.
        """
        return np.array(
            [self._parameters.fsky[i] for i in range(self._parameters.fsky_len)]
        )

    @fsky.setter
    def fsky(self, value):
        try:
            _ = iter(value)
        except TypeError as err:
            raise TypeError(f'The value {value} is not iterable') from err

        try:
            temp = [float(_) for _ in value]
        except TypeError as err:
            raise TypeError('Cannot convert all values to floats') from err

        if not all(_> 0 and _ <= 1 for _ in temp):
            raise ValueError

        if self._parameters.fsky_len:
            free(self._parameters.fsky)

        self._parameters.fsky_len = len(temp)
        self._parameters.fsky = <double *> malloc(sizeof(double) * len(temp))
        for i in range(self._parameters.fsky_len):
            self._parameters.fsky[i] = temp[i]
        self._free_corrfunc()
        self._free_multipoles()
        self._free_covariance_multipoles()


    @property
    def deltaz(self):
        """
        The half width of the redshift bin which is centered at `z_mean`.
        """
        return np.array(
            [self._parameters.deltaz[i] for i in range(self._parameters.deltaz_len)]
        )

    @deltaz.setter
    def deltaz(self, value):
        try:
            _ = iter(value)
        except TypeError as err:
            raise TypeError(f'The value {value} is not iterable') from err

        try:
            temp = [float(_) for _ in value]
        except TypeError as err:
            raise TypeError('Cannot convert all values to floats') from err

        if not all(_>= 0 and _ < 10 for _ in temp):
            raise ValueError

        if self._parameters.deltaz_len:
            free(self._parameters.deltaz)

        self._parameters.deltaz_len = len(temp)
        self._parameters.deltaz = <double *> malloc(sizeof(double) * len(temp))
        for i in range(self._parameters.deltaz_len):
            self._parameters.deltaz[i] = temp[i]
        self._free_corrfunc()
        self._free_multipoles()
        self._free_covariance_multipoles()


    @property
    def integration_sampling(self):
        """
        Returns the maximum number of points to be sampled when using multidimensional
        integration.
        """
        return self._parameters.integration_bins

    @integration_sampling.setter
    def integration_sampling(self, value):
        # the latter is at least the size of a long
        _check_parameter('integration_sampling', value, int, 0, self._max_size)

        self._parameters.integration_bins = value


    @property
    def background_sampling(self):
        """
        Returns the number of points we use to sample the background (from z = 0 to z = 15).
        """
        return self._parameters.background_bins

    @background_sampling.setter
    def background_sampling(self, value):
        _check_parameter('background_sampling', value, int, 0, self._max_size)

        self._parameters.background_bins = value
        self._free_except_parameters()


    @property
    def has_binned_covariance(self):
        r"""
        Returns whether or not the covariance should be bin averaged (see eq.
        (A18) of [arXiv:1509.04293](https://arxiv.org/abs/1509.04293)). Note
        that the parameter `pixelsize` controls the bin width (in
        $\mathrm{Mpc}$) in each redshift bin. If set to False, the covariance
        will not be bin averaged, and eq. (2.52) from
        [arXiv:1806.11090](https://arxiv.org/abs/1806.11090) will be used.
        """
        return bool(self._parameters.covariance_window)

    @has_binned_covariance.setter
    def has_binned_covariance(self, value):
        self._parameters.covariance_window = int(bool(value))
        self._free_covariance_multipoles()


    @property
    def has_density(self):
        """
        Returns whether the density contribution is taken into account.
        """
        return bool(self._parameters.correlation_contrib.den)

    @has_density.setter
    def has_density(self, value : bool):
        self._parameters.correlation_contrib.den = int(bool(value))
        self._free_corrfunc()
        self._free_multipoles()
        self._free_covariance_multipoles()


    @property
    def has_rsd(self):
        """
        Returns whether the RSD contribution is taken into account.
        """
        return bool(self._parameters.correlation_contrib.rsd)

    # TODO implement caching if we're changing to the same value
    @has_rsd.setter
    def has_rsd(self, value : bool):
        self._parameters.correlation_contrib.rsd = int(bool(value))
        self._free_corrfunc()
        self._free_multipoles()
        self._free_covariance_multipoles()


    @property
    def has_lensing(self):
        """
        Returns whether the lensing contribution is taken into account.
        """
        return bool(self._parameters.correlation_contrib.len)

    @has_lensing.setter
    def has_lensing(self, value : bool):
        self._parameters.correlation_contrib.len = int(bool(value))
        self._free_corrfunc()
        self._free_multipoles()


    @property
    def has_d1(self):
        """
        Returns whether the Doppler 1 contribution is taken into account.
        """
        return bool(self._parameters.correlation_contrib.d1)

    @has_d1.setter
    def has_d1(self, value : bool):
        self._parameters.correlation_contrib.d1 = int(bool(value))
        self._free_except_parameters()


    @property
    def has_d2(self):
        """
        Returns whether the Doppler 2 contribution is taken into account.
        """
        return bool(self._parameters.correlation_contrib.d2)

    @has_d2.setter
    def has_d2(self, value : bool):
        self._parameters.correlation_contrib.d2 = int(bool(value))
        self._free_except_parameters()


    @property
    def has_g1(self):
        """
        Returns whether the relativistic non-integrated contribution 1 is taken
        into account.
        """
        return bool(self._parameters.correlation_contrib.g1)

    @has_g1.setter
    def has_g1(self, value : bool):
        self._parameters.correlation_contrib.g1 = int(bool(value))
        self._free_except_parameters()


    @property
    def has_g2(self):
        """
        Returns whether the relativistic non-integrated contribution 2 is taken
        into account.
        """
        return bool(self._parameters.correlation_contrib.g2)

    @has_g2.setter
    def has_g2(self, value : bool):
        self._parameters.correlation_contrib.g2 = int(bool(value))
        self._free_except_parameters()


    @property
    def has_g3(self):
        """
        Returns whether the relativistic non-integrated contribution 3 is taken
        into account.
        """
        return bool(self._parameters.correlation_contrib.g3)

    @has_g3.setter
    def has_g3(self, value : bool):
        self._parameters.correlation_contrib.g3 = int(bool(value))
        self._free_except_parameters()


    @property
    def has_g4(self):
        """
        Returns whether the relativistic integrated contribution 1 is taken
        into account.
        """
        return bool(self._parameters.correlation_contrib.g4)

    @has_g4.setter
    def has_g4(self, value : bool):
        self._parameters.correlation_contrib.g4 = int(bool(value))
        self._free_except_parameters()


    @property
    def has_g5(self):
        """
        Returns whether the relativistic integrated contribution 2 is taken
        into account.
        """
        return bool(self._parameters.correlation_contrib.g5)

    @has_g5.setter
    def has_g5(self, value : bool):
        self._parameters.correlation_contrib.g5 = int(bool(value))
        self._free_except_parameters()


    def reset_contributions(self):
        """
        Helper function that resets all writable `has_*` attributes to False.
        """
        for prop in self.parameters:
            if 'has_' in prop:
                setattr(self, prop, False)


    @property
    def has_flatsky_local(self):
        """
        Whether the flat-sky approximation should be used for local terms (density, RSD).
        """
        return bool(self._parameters.flatsky_local)

    @has_flatsky_local.setter
    def has_flatsky_local(self, value):
        self._parameters.flatsky_local = int(bool(value))
        self._free_integrals()
        self._free_corrfunc()
        self._free_multipoles()


    @property
    def has_flatsky_local_nonlocal(self):
        """
        Whether the flat-sky approximation should be used for cross-terms
        between local and non-local terms (currently density-lensing and RSD-lensing only).
        """
        return bool(self._parameters.flatsky_local_nonlocal)

    @has_flatsky_local_nonlocal.setter
    def has_flatsky_local_nonlocal(self, value):
        self._parameters.flatsky_local_nonlocal = int(bool(value))
        self._free_integrals()
        self._free_corrfunc()
        self._free_multipoles()


    @property
    def has_flatsky_nonlocal(self):
        """
        Whether the flat-sky approximation should be used for integrated
        (non-local) terms (lensing only).
        """
        return bool(self._parameters.flatsky_nonlocal)

    @has_flatsky_nonlocal.setter
    def has_flatsky_nonlocal(self, value):
        self._parameters.flatsky_nonlocal = int(bool(value))
        self._free_integrals()
        self._free_corrfunc()
        self._free_multipoles()


    def maximum_separation(self, z_mean : float, deltaz : float):
        r"""
        Returns the maximum allowed comoving separation (in $\mathrm{Mpc}$) to
        compute the multipoles for a given redshift bin, assuming the current
        cosmology.

        Parameters
        ----------
        z_mean : float
            the mean redshift

        deltaz : float
            the half-width of the redshift bin

        Returns
        -------
        float
        """
        _check_parameter('z_mean', z_mean, (int, float), 0, 15)
        _check_parameter('deltaz', deltaz, (int, float), 0, z_mean)

        if not self._background.flag:
            self._background_init()

        return 2 * (
            ccoffe.coffe_interp_spline(&self._background.comoving_distance, z_mean + deltaz)
            -
            ccoffe.coffe_interp_spline(&self._background.comoving_distance, z_mean)
        )


    @property
    def pk_type(self):
        """
        Returns the currently set power spectrum type (COFFE linear, CLASS linear,
        nonlinear halofit, nonlinear HMcode)
        Possible values: 'linear', 'linear_class', 'halofit', 'hmcode'
        """
        return _allowed_pk_types[self._parameters.pk_type]

    @pk_type.setter
    def pk_type(self, value):
        # we don't care about the case (user is lazy, and so am I)
        if value.lower() not in _allowed_pk_types.values():
            raise ValueError(
                f'The value \'{value}\' is not one of: {list(_allowed_pk_types.values())}'
            )
        if value != self.pk_type:
            self._parameters.pk_type = _allowed_pk_types_inverse[value]
            self._free_power_spectrum()
            self._free_integrals()


    def cross_spectrum(self, k : float, z1 : float, z2 : float, approximation = "geometric"):
        r"""
        Evaluates the matter cross spectrum (in $\mathrm{Mpc}^3$) at some k and z1 and z2.
        Input k must be in $1/\mathrm{Mpc}$.

        Parameters
        ----------
        k : float
            the wavenumber at which to evaluate the cross spectrum

        z1 : float
            the redshift of the first population of tracers

        z2 : float
            the redshift of the second population of tracers

        approximation : str, optional
            which approximation to use to evaluate the cross-spectrum (default: "geometric")

        Returns
        -------
        float
        """
        _check_parameter('k', k, (int, float), self.k_min, self.k_max)
        _check_parameter('z1', z1, (int, float), 0, 15)
        _check_parameter('z2', z2, (int, float), 0, 15)

        if not self._power_spectrum_flag:
            self._power_spectrum_init()

        if not self._background.flag:
            self._background_init()

        return (
            ccoffe.coffe_interp_spline(&self._parameters.power_spectrum, k)
           *ccoffe.coffe_interp_spline(&self._background.D1, z1)
           *ccoffe.coffe_interp_spline(&self._background.D1, z2)
        )


    def power_spectrum(self, k : float, z : float):
        r"""
        Evaluates the matter power spectrum at some k and z.
        Input k must be in $1/\mathrm{Mpc}$.

        Parameters
        ----------
        k : float
            the wavenumber at which to evaluate the spectrum

        z : float
            the mean redshift at which to evaluate the spectrum

        Returns
        -------
        float
        """
        _check_parameter('k', k, (int, float), self.k_min, self.k_max)
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._power_spectrum_flag:
            self._power_spectrum_init()

        if not self._background.flag:
            self._background_init()

        if self._parameters.pk_type == ccoffe.COFFE_PK_LINEAR:
            return (
                ccoffe.coffe_interp_spline(&self._parameters.power_spectrum, k)
               *ccoffe.coffe_interp_spline(&self._background.D1, z)**2
            )
        return ccoffe.coffe_interp_spline2d(&self._parameters.power_spectrum2d, z, k)


    @property
    def k_min(self):
        r"""
        The minimum wavenumber (in $1/\mathrm{Mpc}$) for which the power
        spectrum should be computed.
        """
        return self._parameters.k_min

    @k_min.setter
    def k_min(self, value):
        _check_parameter('k_min', value, (int, float), 0, 1e-3)
        self._parameters.k_min = value
        self._parameters.k_min_norm = value


    @property
    def k_max(self):
        r"""
        The maximum wavenumber (in $1/\mathrm{Mpc}$) for which the power
        spectrum should be computed.
        """
        return self._parameters.k_max

    @k_max.setter
    def k_max(self, value):
        _check_parameter('k_max', value, (int, float), 0.1, 1e3)
        self._parameters.k_max = value
        self._parameters.k_max_norm = value


    def set_power_spectrum_linear(self, k : List[float], pk : List[float], z : float = 0):
        r"""
        Sets the linear matter power spectrum, optionally at some redshift (by
        default it's assumed the input is at z = 0).
        The input $k$ must be in units $1/\mathrm{Mpc}$, and $P(k)$ in units
        $\mathrm{Mpc}^3$.

        Parameters
        ----------
        k : List[float]
            the list of wavenumbers

        pk : List[float]
            the list of values of the linear power spectrum

        z : float, default = 0
            the redshift at which the power spectrum is evaluated

        Returns
        -------
        None
        """
        # value checking
        if not np.allclose(k, np.sort(k)):
            raise ValueError('The input k array must be sorted')

        if len(k) != len(pk):
            raise ValueError('The input arrays have mismatching lengths')

        # TODO error checking for z and rescaling if z > 0

        self._free_power_spectrum()
        self._parameters.k_min = k[0]
        self._parameters.k_max = k[-1]
        self._parameters.k_min_norm = k[0]
        self._parameters.k_max_norm = k[-1]
        size = len(k)

        cdef double *x = <double *>ccoffe.coffe_malloc(sizeof(double) * size)
        cdef double *y = <double *>ccoffe.coffe_malloc(sizeof(double) * size)
        cdef double *x_norm = <double *>ccoffe.coffe_malloc(sizeof(double) * size)
        cdef double *y_norm = <double *>ccoffe.coffe_malloc(sizeof(double) * size)

        for (i, ki), pki in zip(enumerate(k), pk):
            x[i] = ki
            y[i] = pki
            x_norm[i] = ki
            y_norm[i] = pki

        ccoffe.coffe_init_spline(
            &self._parameters.power_spectrum,
            x, y, size, self._parameters.interp_method
        )

        ccoffe.coffe_init_spline(
            &self._parameters.power_spectrum_norm,
            x_norm, y_norm, size, self._parameters.interp_method
        )

        self._power_spectrum_flag = 1

        free(x)
        free(y)
        free(x_norm)
        free(y_norm)

        self._free_integrals()


    def comoving_distance(self, z : float):
        r"""
        Evaluates the comoving distance at some redshift (in units
        $\mathrm{Mpc}$).

        Parameters
        ----------
        z : float
            the redshift

        Returns
        -------
        float
            the result
        """
        # TODO figure out how to dereference an operator
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._background.flag:
            self._background_init()

        return ccoffe.coffe_interp_spline(&self._background.comoving_distance, z)


    def hubble_rate(self, z : float):
        r"""
        Evaluates the Hubble rate ($H(z)$) at some redshift (in units
        $1/\mathrm{Mpc}$).

        Parameters
        ----------
        z : float
            the redshift

        Returns
        -------
        float
            the result
        """
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._background.flag:
            self._background_init()

        return ccoffe.coffe_interp_spline(&self._background.Hz, z)


    def scale_factor(self, z : float):
        """
        Returns the scale factor $a(z)$ evaluated at some redshift.

        Parameters
        ----------
        z : float
            the redshift

        Returns
        -------
        float
            the result
        """
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._background.flag:
            self._background_init()

        return ccoffe.coffe_interp_spline(&self._background.a, z)


    def hubble_rate_conformal(self, z : float):
        r"""
        Evaluates the conformal Hubble rate ($\mathcal{H}(z)$) at some redshift
        (in units $1/\mathrm{Mpc}$).

        Parameters
        ----------
        z : float
            the redshift

        Returns
        -------
        float
            the result
        """
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._background.flag:
            self._background_init()

        return ccoffe.coffe_interp_spline(&self._background.conformal_Hz, z)


    def hubble_rate_conformal_derivative(self, z : float):
        r"""
        Evaluates the first derivative of the conformal Hubble rate w.r.t.
        conformal time ($\mathrm{d} \mathcal{H}(z)/\mathrm{d}\tau$) at some
        redshift (in units $1/\mathrm{Mpc}^2$).

        Parameters
        ----------
        z : float
            the redshift

        Returns
        -------
        float
            the result
        """
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._background.flag:
            self._background_init()

        return ccoffe.coffe_interp_spline(&self._background.conformal_Hz_prime, z)


    def growth_factor(self, z : float):
        r"""
        Returns the scale-independent function $D_1$ evaluated at some redshift.

        Parameters
        ----------
        z : float
            the redshift

        Returns
        -------
        float
            the result
        """
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._background.flag:
            self._background_init()

        return ccoffe.coffe_interp_spline(&self._background.D1, z)


    def growth_rate(self, z : float):
        r"""
        Returns the scale-independent growth rate $f \equiv \mathrm{d}(\log
        D_1) / \mathrm{d} (\log a)$ evaluated at some redshift.

        Parameters
        ----------
        z : float
            the redshift

        Returns
        -------
        float
            the result
        """
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._background.flag:
            self._background_init()

        return ccoffe.coffe_interp_spline(&self._background.f, z)


    def compute_corrfunc(
        self, *,
        # TODO how should we name the parameters?
        z : float, r : float, mu : float,
        recompute : bool = False,
    ):
        r"""
        Computes the correlation function of the current configuration at the
        point $(\bar{z}, r, \mu)$.

        Parameters
        ----------
        z : float
            the mean redshift

        r : float
            the comoving separation (in $\mathrm{Mpc}$) between the two points
            in the sky

        mu : float
            the angle mu (see `help(coffe.mu)` for definition)

        Returns
        -------
        float
            the result
        """
        self._check_coords_corrfunc()
        self._check_contributions()

        if not self._background.flag or recompute:
            self._background_init()

        if not self._power_spectrum_flag or recompute:
            self._power_spectrum_init()

        if not self._integral.size or recompute:
            self._integrals_init()

        _check_parameter('z', z, (int, float), 0, 15)
        _check_parameter('r', r, (int, float), 0, 1500)
        _check_parameter('mu', mu, (int, float), -1, 1)

        cdef ccoffe.gsl_error_handler_t *default_handler = \
            ccoffe.gsl_set_error_handler_off()

        result = ccoffe.coffe_integrate(
            &self._parameters,
            &self._background,
            &self._integral,
            z, r, mu, 0,
            ccoffe.NONINTEGRATED, ccoffe.CORRFUNC
        ) + \
        ccoffe.coffe_integrate(
            &self._parameters,
            &self._background,
            &self._integral,
            z, r, mu, 0,
            ccoffe.SINGLE_INTEGRATED, ccoffe.CORRFUNC
        ) + \
        ccoffe.coffe_integrate(
            &self._parameters,
            &self._background,
            &self._integral,
            z, r, mu, 0,
            ccoffe.DOUBLE_INTEGRATED, ccoffe.CORRFUNC
        )

        ccoffe.gsl_set_error_handler(default_handler)

        return result


    def compute_multipole(
        self, *,
        z : float, r : float, l : int,
        recompute : bool = False,
    ):
        r"""
        Computes the multipole of the current configuration at the point
        $(\bar{z}, r, \ell)$.

        Parameters
        ----------
        z : float
            the mean redshift

        r : float
            the comoving separation (in $\mathrm{Mpc}$) between the 2 points in
            the sky

        l : int
            the multipole moment

        Returns
        -------
        float
            the result
        """
        self._check_coords_multipoles()
        self._check_contributions()
        if not self._background.flag or recompute:
            self._background_init()

        if not self._power_spectrum_flag or recompute:
            self._power_spectrum_init()

        if not self._integral.size or recompute:
            self._integrals_init()

        _check_parameter('z', z, (int, float), 0, 15)
        _check_parameter('r', r, (int, float), 0, 1500)
        _check_parameter('l', l, (int,), 0, 10)

        cdef ccoffe.gsl_error_handler_t *default_handler = \
            ccoffe.gsl_set_error_handler_off()

        result = ccoffe.coffe_integrate(
            &self._parameters,
            &self._background,
            &self._integral,
            z, r, 0, l,
            ccoffe.NONINTEGRATED, ccoffe.MULTIPOLES
        ) + \
        ccoffe.coffe_integrate(
            &self._parameters,
            &self._background,
            &self._integral,
            z, r, 0, l,
            ccoffe.SINGLE_INTEGRATED, ccoffe.MULTIPOLES
        ) + \
        ccoffe.coffe_integrate(
            &self._parameters,
            &self._background,
            &self._integral,
            z, r, 0, l,
            ccoffe.DOUBLE_INTEGRATED, ccoffe.MULTIPOLES
        )

        ccoffe.gsl_set_error_handler(default_handler)

        return result


    def compute_multipoles_bulk(self, recompute : bool = False):
        """
        Computes the multipoles in bulk for all currently set separations,
        redshifts, and multipole moments.

        Parameters
        ----------
        recompute : bool, default = False
            if set to True, recomputes all of the other modules first to make
            sure we haven't forgotten to refresh one of them. Mainly useful for debugging
            purposes.

        Returns
        -------
        an array of instances of `coffe.representation.Multipoles`.
        """
        self._check_coords_multipoles()
        self._check_contributions()
        recompute = bool(recompute)
        if not self._background.flag or recompute:
            self._background_init()

        if not self._power_spectrum_flag or recompute:
            self._power_spectrum_init()

        if not self._integral.size or recompute:
            self._integrals_init()

        if not self._multipoles.size or recompute:
            self._multipoles_init()

        return np.array([
            Multipoles(
                z=self._multipoles.array[i].coords.z_mean,
                r=self._multipoles.array[i].coords.separation,
                l=self._multipoles.array[i].coords.l,
                value=self._multipoles.array[i].value,
            ) for i in range(self._multipoles.size)
        ])


    def compute_corrfunc_bulk(self, recompute : bool = False):
        """
        Computes the 2PCF in bulk for all currently set separations,
        redshifts, and angles.

        Parameters
        ----------
        recompute : bool, default = False
            if set to True, recomputes all of the other modules first to make
            sure we haven't forgotten to refresh one of them. Mainly useful for debugging
            purposes.

        Returns
        -------
        an array of instances of `coffe.representation.Corrfunc`.

        """
        self._check_coords_corrfunc()
        self._check_contributions()
        recompute = bool(recompute)
        if not self._background.flag or recompute:
            self._background_init()

        if not self._power_spectrum_flag or recompute:
            self._power_spectrum_init()

        if not self._integral.size or recompute:
            self._integrals_init()

        if not self._corrfunc.size or recompute:
            self._corrfunc_init()

        return np.array([
            Corrfunc(
                z=self._corrfunc.array[i].coords.z_mean,
                r=self._corrfunc.array[i].coords.separation,
                mu=self._corrfunc.array[i].coords.mu,
                value=self._corrfunc.array[i].value,
            ) for i in range(self._corrfunc.size)
        ])


    def compute_covariance_bulk(self, recompute : bool = False):
        """
        Computes the covariance of the multipoles in bulk for all currently set
        separations, redshifts, and multipole moments.

        Parameters
        ----------
        recompute : bool, default = False
            if set to True, recomputes all of the other modules first to make
            sure we haven't forgotten to refresh one of them. Mainly useful for debugging
            purposes.

        Returns
        -------
        an array of instances of `coffe.representation.Covariance`.
        """
        self._check_coords_multipoles()
        self._check_contributions()
        recompute = bool(recompute)

        if not all(
            len(self.z_mean) == len(_)
            for _ in (self.deltaz, self.number_density1, self.number_density2, self.pixelsize, self.fsky)
        ):
            raise ValueError('Mismatching lengths for covariance parameters (`z_mean`, `deltaz`, `number_density1`, `number_density2`, `pixelsize`, `fsky`')

        if not self._power_spectrum_flag or recompute:
            self._power_spectrum_init()

        if not self._background.flag or recompute:
            self._background_init()

        try:
            [self.maximum_separation(zi, deltazi) for zi, deltazi in zip(self.z_mean, self.deltaz)]
        except ValueError as err:
            raise ValueError from err

        temp = self._parameters.output_type
        self._parameters.output_type = ccoffe.COVARIANCE_MULTIPOLES
        if not self._covariance_multipoles.size or recompute:
            self._covariance_init()
        self._parameters.output_type = temp

        return np.array([
            Covariance(
                z=self._covariance_multipoles.array[i].coords.z_mean,
                r1=self._covariance_multipoles.array[i].coords.separation1,
                r2=self._covariance_multipoles.array[i].coords.separation2,
                l1=self._covariance_multipoles.array[i].coords.l1,
                l2=self._covariance_multipoles.array[i].coords.l2,
                value=self._covariance_multipoles.array[i].value,
            ) for i in range(self._covariance_multipoles.size)
        ])


    def _background_init(self):
        ccoffe.coffe_background_init(
            &self._parameters,
            &self._background
        )

    def _power_spectrum_init(self):
        ccoffe.parse_external_power_spectrum(&self._parameters)
        self._power_spectrum_flag = 1

    def _integrals_init(self):
        ccoffe.coffe_integrals_init(
            &self._parameters,
            &self._background,
            &self._integral
        )

    def _corrfunc_init(self):
        ccoffe.coffe_corrfunc_init(
            &self._parameters,
            &self._background,
            &self._integral,
            &self._corrfunc
        )

    def _multipoles_init(self):
        ccoffe.coffe_multipoles_init(
            &self._parameters,
            &self._background,
            &self._integral,
            &self._multipoles
        )

    def _covariance_init(self):
        ccoffe.coffe_covariance_init(
            &self._parameters,
            &self._background,
            &self._covariance_multipoles,
            &self.__dummy
        )
