# distutils: sources = src/errors.c src/common.c src/parser.c src/background.c src/twofast.c src/integrals.c src/signal.c src/functions.c src/corrfunc.c src/multipoles.c src/utils.c src/twobessel.c src/covariance.c
# distutils: include_dirs = src/ ./
# distutils: libraries = m gsl gslcblas config fftw3 cuba class gomp
# distutils: extra_compile_args = ['-fopenmp', '-Ofast']
# distutils: extra_link_args = ['-fopenmp']

# TODO figure out how to use OpenMP

from libc.stdlib cimport malloc, free
cimport ccoffe

_COFFE_HUBBLE = (1./(2997.92458))


from cython.operator cimport dereference
from ctypes import CFUNCTYPE

from typing import Any, Callable, List, Tuple, Union
from dataclasses import dataclass
import numpy as np


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
        raise TypeError
    if xmin is not None and value < xmin:
            raise ValueError
    if xmax is not None and value > xmax:
            raise ValueError



class Representation:
    def to_dict(self):
        """
        The representation of the class as a dictionary.
        """
        return {
            key : getattr(self, key) \
            for key in dir(self.__class__) \
            if hasattr(getattr(self.__class__, key), '__set__') \
            and not key.startswith('__')
        }


    def __repr__(self):
        """
        User-friendly representation of the class
        """
        return f'{self.__class__}({self.to_dict()})'


    def _repr_html_(self):
        names = self.to_dict().keys()
        values = self.to_dict().values()
        temp = (
            '<tr>' + ('<th>{}</th>' * len(names)).format(*names) + '</tr>'
        ) if names else ''
        header = f'<thead>{temp}</thead>'

        body = '<tbody>' + (
            '<td>{}</td>' * len(values)
        ).format(*values) + '</tbody>'

        return f'<table>{header}{body}</table>'



class Covariance(Representation):
    def __init__(
        self, *,
        r1 : float, r2 : float,
        l1 : int, l2 : int,
        z : float,
        value : float,
    ):
        self._r1 = r1
        self._r2 = r2
        self._l1 = l1
        self._l2 = l2
        self._z = z
        self._value = value

    @property
    def r1(self):
        return self._r1

    @property
    def r2(self):
        return self._r2

    @property
    def l1(self):
        return self._l1

    @property
    def l2(self):
        return self._l2

    @property
    def z(self):
        return self._z

    @property
    def value(self):
        return self._value



class Corrfunc(Representation):
    def __init__(
        self, *,
        r : float, mu : float, z : float,
        value : float,
    ):
        self._r = r
        self._mu = mu
        self._z = z
        self._value = value

    @property
    def mu(self):
        return self._mu

    @property
    def r(self):
        return self._r

    @property
    def z(self):
        return self._z

    @property
    def value(self):
        return self._value



class Multipoles(Representation):
    def __init__(
        self, *,
        l : int, r : float, z : float,
        value : float,
    ):
        self._l = l
        self._r = r
        self._z = z
        self._value = value

    @property
    def l(self):
        return self._l

    @property
    def r(self):
        return self._r

    @property
    def z(self):
        return self._z

    @property
    def value(self):
        return self._value



cdef class Coffe:
    cdef ccoffe.coffe_parameters_t _parameters
    cdef ccoffe.coffe_background_t _background
    cdef ccoffe.coffe_integral_array_t _integral
    cdef ccoffe.coffe_corrfunc_array_t _corrfunc
    cdef ccoffe.coffe_multipoles_array_t _multipoles
    cdef ccoffe.coffe_covariance_array_t _covariance_multipoles
    cdef ccoffe.coffe_covariance_array_t __dummy
    cdef int _power_spectrum_flag

    # TODO make it possible to set parameters directly in the constructor
    def __cinit__(self):
        """
        Constructor that initializes the structures and sets default parameters.
        """
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


    def _free_background(self):
        ccoffe.coffe_background_free(&self._background)


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
        self._free_corrfunc()
        self._free_multipoles()
        self._free_covariance_multipoles()


    def __dealloc__(self):
        self._free_except_parameters()
        ccoffe.coffe_parameters_free(&self._parameters)


    def _balance_content(self):
        """
        Balances the energy content of dark energy so it all adds up to 1
        """
        self._parameters.Omega0_de = 1 - self._parameters.Omega0_cdm - \
            self._parameters.Omega0_baryon - self._parameters.Omega0_gamma


    def parameters(self):
        """
        Returns the current writable parameters as a dictionary (can be
        re-used for `set_parameters`).
        """
        properties = {
            key : getattr(self, key) \
            for key in dir(self.__class__) \
            if hasattr(getattr(self.__class__, key), '__set__')
        }
        writable = {}
        for key in properties:
            try:
                setattr(self, key, getattr(self, key))
                writable[key] = getattr(self, key)
            except AttributeError:
                pass
        return writable


    def set_parameters(self, value : dict):
        """
        Bulk setter of parameters.
        The passed object _must_ be a dictionary with strings as keys.
        """
        if not isinstance(value, dict):
            raise TypeError

        for key in value:
            if not hasattr(self, key):
                raise AttributeError(
                    f'The parameter {key} is not a valid parameter in COFFE.'
                )

        for key in value:
            setattr(self, key, value[key])


    # TODO how do we update the other omegas here (like omega_m)?
    @property
    def omega_cdm(self):
        return self._parameters.Omega0_cdm

    @omega_cdm.setter
    def omega_cdm(self, value):
        if not isinstance(value, (int, float)):
            raise TypeError
        if value <= 0 or value >= 1:
            raise ValueError
        if not np.allclose(value, self.omega_cdm):
            # we set the value, rebalance the Omega budget, and free memory
            self._parameters.Omega0_cdm = value
            self._balance_content()
            self._free_except_parameters()


    @property
    def omega_baryon(self):
        return self._parameters.Omega0_baryon

    @omega_baryon.setter
    def omega_baryon(self, value):
        if not isinstance(value, (int, float)):
            raise TypeError
        if value <= 0 or value >= 1:
            raise ValueError
        if not np.allclose(value, self.omega_baryon):
            # we set the value, rebalance the Omega budget, and free memory
            self._parameters.Omega0_baryon = value
            self._balance_content()
            self._free_except_parameters()


    @property
    def omega_gamma(self):
        return self._parameters.Omega0_gamma

    @omega_gamma.setter
    def omega_gamma(self, value):
        _check_parameter('omega_gamma', value, (int, float), 0, 1)
        if not np.allclose(value, self.omega_gamma):
            self._parameters.Omega0_gamma = value
            self._balance_content()
            self._free_except_parameters()


    @property
    def omega_de(self):
        return self._parameters.Omega0_de


    @property
    def h(self):
        return self._parameters.h

    @h.setter
    def h(self, value):
        if not isinstance(value, (int, float)):
            raise TypeError
        if value <= 0 or value >= 1:
            raise ValueError
        if not np.allclose(value, self.h):
            # we set the value, but don't free the background since it's unaffected by h
            self._parameters.h = value
            self._free_integrals()
            self._free_corrfunc()
            self._free_multipoles()
            self._free_covariance_multipoles()


    @property
    def w0(self):
        return self._parameters.w0

    @w0.setter
    def w0(self, value):
        _check_parameter('w0', value, (int, float), -2, 0)
        if not np.allclose(value, self.w0):
            self._parameters.w0 = value
            self._free_except_parameters()


    @property
    def wa(self):
        return self._parameters.wa

    @wa.setter
    def wa(self, value):
        _check_parameter('wa', value, (int, float), -1, 1)
        if not np.allclose(value, self.wa):
            self._parameters.wa = value
            self._free_except_parameters()


    @property
    def n_s(self):
        return self._parameters.n_s

    @n_s.setter
    def n_s(self, value):
        _check_parameter('n_s', value, (int, float), 0.5, 1.5)
        if not np.allclose(value, self.n_s):
            self._parameters.n_s = value
            self._free_integrals()
            self._free_corrfunc()
            self._free_multipoles()
            self._free_covariance_multipoles()


    @property
    def sigma8(self):
        return self._parameters.sigma8

    @sigma8.setter
    def sigma8(self, value):
        _check_parameter('sigma8', value, (int, float), 0, 2)
        if not np.allclose(value, self.sigma8):
            self._parameters.sigma8 = value
            self._free_integrals()
            self._free_corrfunc()
            self._free_multipoles()
            self._free_covariance_multipoles()


    @property
    def sep(self):
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
        return np.array(
            [self._parameters.mu[i] for i in range(self._parameters.mu_len)]
        )


    def galaxy_bias1(self, z : float):
        """
        Evaluates the galaxy bias of the first population at some redshift.
        """
        # TODO figure out how to dereference an operator
        _check_parameter('z', z, (int, float), 0, 15)

        return ccoffe.coffe_interp_spline(&self._parameters.galaxy_bias1, z)


    def set_galaxy_bias1(self, value : Callable):
        """
        Sets the value of the galaxy bias.
        The value set must be callable.
        """
        x_sampling = np.linspace(0, 15, 1000)
        cdef double *x = NULL
        cdef double *y = NULL

        try:
            [value(_) for _ in x_sampling]
        except TypeError as err:
            raise TypeError(f'Unable to sample function {value}') from err

        x = <double *>malloc(sizeof(double) * len(x_sampling))
        y = <double *>malloc(sizeof(double) * len(x_sampling))
        for i in range(len(x_sampling)):
            x[i] = x_sampling[i]
            y[i] = value(x[i])
        ccoffe.coffe_init_spline(
            &self._parameters.galaxy_bias1,
            x, y, len(x_sampling), self._parameters.interp_method
        )

        if (x != NULL):
            free(x)
            x = NULL
        if (y != NULL):
            free(y)
            y = NULL

        self._free_corrfunc()
        self._free_multipoles()


    def galaxy_bias2(self, z : float):
        """
        Evaluates the galaxy bias of the first population at some redshift.
        """
        # TODO figure out how to dereference an operator
        _check_parameter('z', z, (int, float), 0, 15)

        return ccoffe.coffe_interp_spline(&self._parameters.galaxy_bias2, z)


    def set_galaxy_bias2(self, value : Callable):
        """
        Sets the value of the galaxy bias.
        The value set must be callable.
        """
        x_sampling = np.linspace(0, 15, 1000)
        cdef double *x = NULL
        cdef double *y = NULL

        try:
            [value(_) for _ in x_sampling]
        except TypeError as err:
            raise TypeError(f'Unable to sample function {value}') from err

        x = <double *>malloc(sizeof(double) * len(x_sampling))
        y = <double *>malloc(sizeof(double) * len(x_sampling))
        for i in range(len(x_sampling)):
            x[i] = x_sampling[i]
            y[i] = value(x[i])
        ccoffe.coffe_init_spline(
            &self._parameters.galaxy_bias2,
            x, y, len(x_sampling), self._parameters.interp_method
        )

        if (x != NULL):
            free(x)
            x = NULL
        if (y != NULL):
            free(y)
            y = NULL

        self._free_corrfunc()
        self._free_multipoles()


    def magnification_bias1(self, z : float):
        """
        Evaluates the galaxy bias of the first population at some redshift.
        """
        # TODO figure out how to dereference an operator
        _check_parameter('z', z, (int, float), 0, 15)

        return ccoffe.coffe_interp_spline(&self._parameters.magnification_bias1, z)


    def set_magnification_bias1(self, value : Callable):
        """
        Sets the value of the galaxy bias.
        The value set must be callable.
        """
        x_sampling = np.linspace(0, 15, 1000)
        cdef double *x = NULL
        cdef double *y = NULL

        try:
            [value(_) for _ in x_sampling]
        except TypeError as err:
            raise TypeError(f'Unable to sample function {value}') from err

        x = <double *>malloc(sizeof(double) * len(x_sampling))
        y = <double *>malloc(sizeof(double) * len(x_sampling))
        for i in range(len(x_sampling)):
            x[i] = x_sampling[i]
            y[i] = value(x[i])
        ccoffe.coffe_init_spline(
            &self._parameters.magnification_bias1,
            x, y, len(x_sampling), self._parameters.interp_method
        )

        if (x != NULL):
            free(x)
            x = NULL
        if (y != NULL):
            free(y)
            y = NULL

        self._free_corrfunc()
        self._free_multipoles()


    def magnification_bias2(self, z : float):
        """
        Evaluates the galaxy bias of the first population at some redshift.
        """
        # TODO figure out how to dereference an operator
        _check_parameter('z', z, (int, float), 0, 15)

        return ccoffe.coffe_interp_spline(&self._parameters.magnification_bias2, z)


    def set_magnification_bias2(self, value : Callable):
        """
        Sets the value of the galaxy bias.
        The value set must be callable.
        """
        x_sampling = np.linspace(0, 15, 1000)
        cdef double *x = NULL
        cdef double *y = NULL

        try:
            [value(_) for _ in x_sampling]
        except TypeError as err:
            raise TypeError(f'Unable to sample function {value}') from err

        x = <double *>malloc(sizeof(double) * len(x_sampling))
        y = <double *>malloc(sizeof(double) * len(x_sampling))
        for i in range(len(x_sampling)):
            x[i] = x_sampling[i]
            y[i] = value(x[i])
        ccoffe.coffe_init_spline(
            &self._parameters.magnification_bias2,
            x, y, len(x_sampling), self._parameters.interp_method
        )

        if (x != NULL):
            free(x)
            x = NULL
        if (y != NULL):
            free(y)
            y = NULL

        self._free_corrfunc()
        self._free_multipoles()


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
            free(self._parameters.sep)

        self._parameters.mu_len = len(temp)
        self._parameters.mu = <double *> malloc(sizeof(double) * len(temp))
        for i in range(self._parameters.mu_len):
            self._parameters.mu[i] = temp[i]
        # we need to re-compute it since we changed the values
        # TODO make it so that if we give identical values, we don't recompute it
        self._free_corrfunc()


    @property
    def l(self):
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
        self._free_multipoles()
        self._free_covariance_multipoles()


    @property
    def z_mean(self):
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
    def has_flatsky_local(self):
        return bool(self._parameters.flatsky_local)

    @has_flatsky_local.setter
    def has_flatsky_local(self, value):
        self._parameters.flatsky_local = int(bool(value))
        self._free_corrfunc()
        self._free_multipoles()


    @property
    def has_flatsky_local_nonlocal(self):
        return bool(self._parameters.flatsky_local_nonlocal)

    @has_flatsky_local_nonlocal.setter
    def has_flatsky_local_nonlocal(self, value):
        self._parameters.flatsky_local_nonlocal = int(bool(value))
        self._free_corrfunc()
        self._free_multipoles()


    @property
    def has_flatsky_nonlocal(self):
        return bool(self._parameters.flatsky_nonlocal)

    @has_flatsky_nonlocal.setter
    def has_flatsky_nonlocal(self, value):
        self._parameters.flatsky_nonlocal = int(bool(value))
        self._free_corrfunc()
        self._free_multipoles()


    def maximum_separation(self, z_mean : float, deltaz : float):
        """
        Returns the maximum allowed comoving separation (in Mpc/h) to compute the
        multipoles for a given redshift bin, assuming the current cosmology.
        """
        _check_parameter('z_mean', z_mean, (int, float), 0, 15)
        _check_parameter('deltaz', deltaz, (int, float), 0, z_mean)

        if not self._background.flag:
            self._background_init()

        return 2 * (
            ccoffe.coffe_interp_spline(&self._background.comoving_distance, z_mean + deltaz) \
            - \
            ccoffe.coffe_interp_spline(&self._background.comoving_distance, z_mean)
        ) / _COFFE_HUBBLE


    def _power_spectrum_init(self):
        ccoffe.parse_external_power_spectrum(&self._parameters)
        self._power_spectrum_flag = 1


    def cross_spectrum(self, k : float, z1 : float, z2 : float):
        """
        Evaluates the (for now linear) matter cross spectrum at some k and z1 and z2.
        """
        _check_parameter('k', k, (int, float), 1e-5, 1e3)
        _check_parameter('z1', z1, (int, float), 0, 15)
        _check_parameter('z2', z2, (int, float), 0, 15)

        if not self._power_spectrum_flag:
            self._power_spectrum_init()

        if not self._background.flag:
            self._background_init()

        return \
            ccoffe.coffe_interp_spline(&self._parameters.power_spectrum, k) \
           *ccoffe.coffe_interp_spline(&self._background.D1, z1) \
           *ccoffe.coffe_interp_spline(&self._background.D1, z2)


    def power_spectrum(self, k : float, z : float):
        """
        Evaluates the (for now linear) matter power spectrum at some k and z.
        """
        _check_parameter('k', k, (int, float), 1e-5, 1e3)
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._power_spectrum_flag:
            self._power_spectrum_init()

        if not self._background.flag:
            self._background_init()

        return \
            ccoffe.coffe_interp_spline(&self._parameters.power_spectrum, k) \
           *ccoffe.coffe_interp_spline(&self._background.D1, z) \
           *ccoffe.coffe_interp_spline(&self._background.D1, z)


    def comoving_distance(self, z : float):
        """
        Evaluates the comoving distance at some redshift.
        """
        # TODO figure out how to dereference an operator
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._background.flag:
            self._background_init()

        return ccoffe.coffe_interp_spline(&self._background.comoving_distance, z)


    def scale_factor(self, z : float):
        """
        Returns the scale factor a evaluated at some redshift.
        """
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._background.flag:
            self._background_init()

        return ccoffe.coffe_interp_spline(&self._background.a, z)


    def growth_factor(self, z : float):
        """
        Returns the scale-independent function D1 evaluated at some redshift.
        """
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._background.flag:
            self._background_init()

        return ccoffe.coffe_interp_spline(&self._background.D1, z)


    def growth_rate(self, z : float):
        """
        Returns the scale-independent growth rate f evaluated at some redshift.
        """
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._background.flag:
            self._background_init()

        return ccoffe.coffe_interp_spline(&self._background.f, z)


    def compute_corrfunc(
        self, *,
        # TODO how should we name the parameters?
        z : float, r : float, mu : float,
    ):
        """
        Computes the correlation function of the current configuration at the point (z, r, mu)
        """
        if not self._background.flag:
            self._background_init()

        if not self._power_spectrum_flag:
            self._power_spectrum_init()

        if not self._integral.size:
            self._integrals_init()

        _check_parameter('z', z, (int, float), 0, 15)
        _check_parameter('r', r, (int, float), 0, 1500)
        _check_parameter('mu', mu, (int, float), -1, 1)

        return ccoffe.coffe_integrate(
            &self._parameters,
            &self._background,
            &self._integral,
            z, r * _COFFE_HUBBLE, mu, 0,
            ccoffe.NONINTEGRATED, ccoffe.CORRFUNC
        ) + \
        ccoffe.coffe_integrate(
            &self._parameters,
            &self._background,
            &self._integral,
            z, r * _COFFE_HUBBLE, mu, 0,
            ccoffe.SINGLE_INTEGRATED, ccoffe.CORRFUNC
        ) + \
        ccoffe.coffe_integrate(
            &self._parameters,
            &self._background,
            &self._integral,
            z, r * _COFFE_HUBBLE, mu, 0,
            ccoffe.DOUBLE_INTEGRATED, ccoffe.CORRFUNC
        )


    def compute_multipole(
        self, *,
        z : float, r : float, l : int,
    ):
        """
        Computes the multipole of the current configuration at the point (z, r, l)
        """
        if not self._background.flag:
            self._background_init()

        if not self._power_spectrum_flag:
            self._power_spectrum_init()

        if not self._integral.size:
            self._integrals_init()

        _check_parameter('z', z, (int, float), 0, 15)
        _check_parameter('r', r, (int, float), 0, 1500)
        _check_parameter('l', l, (int,), 0, 10)

        return ccoffe.coffe_integrate(
            &self._parameters,
            &self._background,
            &self._integral,
            z, r * _COFFE_HUBBLE, 0, l,
            ccoffe.NONINTEGRATED, ccoffe.MULTIPOLES
        ) + \
        ccoffe.coffe_integrate(
            &self._parameters,
            &self._background,
            &self._integral,
            z, r * _COFFE_HUBBLE, 0, l,
            ccoffe.SINGLE_INTEGRATED, ccoffe.MULTIPOLES
        ) + \
        ccoffe.coffe_integrate(
            &self._parameters,
            &self._background,
            &self._integral,
            z, r * _COFFE_HUBBLE, 0, l,
            ccoffe.DOUBLE_INTEGRATED, ccoffe.MULTIPOLES
        )


    def compute_multipoles_bulk(self):
        """
        Returns whatever `coffe_multipoles_init` returns.
        """
        if not self._background.flag:
            self._background_init()

        if not self._power_spectrum_flag:
            self._power_spectrum_init()

        if not self._integral.size:
            self._integrals_init()

        if not self._multipoles.size:
            self._multipoles_init()

        return np.array([
            Multipoles(
                z=self._multipoles.array[i].coords.z_mean,
                r=self._multipoles.array[i].coords.separation,
                l=self._multipoles.array[i].coords.l,
                value=self._multipoles.array[i].value,
            ) for i in range(self._multipoles.size)
        ])


    def compute_corrfunc_bulk(self):
        """
        Returns whatever `coffe_corrfunc_init` returns.
        """
        if not self._background.flag:
            self._background_init()

        if not self._power_spectrum_flag:
            self._power_spectrum_init()

        if not self._integral.size:
            self._integrals_init()

        if not self._corrfunc.size:
            self._corrfunc_init()

        return np.array([
            Corrfunc(
                z=self._corrfunc.array[i].coords.z_mean,
                r=self._corrfunc.array[i].coords.separation,
                mu=self._corrfunc.array[i].coords.mu,
                value=self._corrfunc.array[i].value,
            ) for i in range(self._corrfunc.size)
        ])


    def compute_covariance_bulk(self):
        """
        Returns whatever `coffe_covariance_init` returns.
        """
        if not self._power_spectrum_flag:
            self._power_spectrum_init()

        if not self._background.flag:
            self._background_init()

        temp = self._parameters.output_type
        self._parameters.output_type = 4
        if not self._covariance_multipoles.size:
            self._covariance_init()
        self._parameters.output_type = temp

        return np.array([
            Covariance(
                z=self._covariance_multipoles.array[i].coords.z_mean,
                r1=self._covariance_multipoles.array[i].coords.separation1,
                r2=self._covariance_multipoles.array[i].coords.separation2,
                l1=self._covariance_multipoles.array[i].coords.l1,
                l2=self._covariance_multipoles.array[i].coords.l2,
                mu=self._covariance_multipoles.array[i].coords.l,
                value=self._covariance_multipoles.array[i].value,
            ) for i in range(self._covariance_multipoles.size)
        ])


    def _background_init(self):
        ccoffe.coffe_background_init(
            &self._parameters,
            &self._background
        )

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
