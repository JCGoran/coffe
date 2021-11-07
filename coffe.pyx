# distutils: sources = src/errors.c src/common.c src/parser.c src/background.c src/twofast.c src/integrals.c src/signal.c src/functions.c src/corrfunc.c src/multipoles.c src/utils.c src/twobessel.c src/covariance.c
# distutils: include_dirs = src/ ./
# distutils: libraries = m gsl gslcblas config fftw3 cuba class
# distutils: extra_compile_flags += ['-fopenmp', '-O3']

# TODO figure out how to use OpenMP

cimport ccoffe

_COFFE_HUBBLE = (1./(2997.92458))


from cython.operator cimport dereference

from typing import Any, List, Tuple, Union
from dataclasses import dataclass
import numpy as np


def _check_parameter(
    name : str,
    value : Any,
    kind : Union[type, Tuple[type], List[type]],
    xmin : float = None,
    xmax : float = None
):
    if not isinstance(value, kind):
        raise TypeError
    if xmin:
        if value < xmin:
            raise ValueError
    if xmax:
        if value > xmax:
            raise ValueError


class Covariance:
    def __init__(
        self, *,
        r1 : float,
        r2 : float,
        l1 : int,
        l2 : int,
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

    def to_dict(self):
        return {
            'r1' : self.r1,
            'r2' : self.r2,
            'l1' : self.l1,
            'l2' : self.l2,
            'z' : self.z,
            'value' : self.value,
        }


class Corrfunc:
    def __init__(
        self, *,
        r : float,
        mu : float,
        z : float,
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

    def to_dict(self):
        return {
            'mu' : self.mu,
            'r' : self.r,
            'z' : self.z,
            'value' : self.value,
        }


class Multipoles:
    def __init__(
        self, *,
        l : int,
        r : float,
        z : float,
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

    def to_dict(self):
        return {
            'l' : self.l,
            'r' : self.r,
            'z' : self.z,
            'value' : self.value,
        }


class Parameter:
    def __init__(
        self,
        kind : type,
        # whether COFFE treats this parameter as a pointer or not (i.e. any kind of array)
        pointer : bool = False,
        xmin : Union[float, None] = None,
        xmax : Union[float, None] = None,
    ):
        pass



cdef class Coffe:
    cdef ccoffe.coffe_parameters_t _parameters
    cdef ccoffe.coffe_background_t _background
    cdef ccoffe.coffe_integral_array_t _integral
    cdef ccoffe.coffe_corrfunc_array_t _corrfunc
    cdef ccoffe.coffe_multipoles_array_t _multipoles
    cdef ccoffe.coffe_covariance_array_t _covariance_multipoles
    cdef ccoffe.coffe_covariance_array_t __dummy

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

        # flag which keeps track of whether a given setting has changed
        # TODO maybe we should have flags which only re-run certain functions;
        # for instance, if we change separations, in principle we don't need to
        # re-run the background module

    def __dealloc__(self):
        ccoffe.coffe_parameters_free(&self._parameters)
        ccoffe.coffe_background_free(&self._background)
        ccoffe.coffe_corrfunc_free(&self._corrfunc)
        ccoffe.coffe_multipoles_free(&self._multipoles)

    def _balance_content(self):
        """
        Balances the energy content of dark energy so it all adds up to 1
        """
        self._parameters.Omega0_de = 1 - self._parameters.Omega0_cdm - \
            self._parameters.Omega0_baryon - self._parameters.Omega0_gamma

    def _background_init(self):
        """
        Runs the background module.
        """
        if self._background.flag:
            ccoffe.coffe_background_free(&self._background)
        ccoffe.coffe_background_init(
            &self._parameters,
            &self._background
        )

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

    @property
    def omega_cdm(self):
        return self._parameters.Omega0_cdm

    @omega_cdm.setter
    def omega_cdm(self, value):
        if not isinstance(value, (int, float)):
            raise TypeError
        if value <= 0 or value >= 1:
            raise ValueError
        self._parameters.Omega0_cdm = value
        self._balance_content()

    @property
    def omega_de(self):
        return self._parameters.Omega0_de

    @property
    def w0(self):
        return self._parameters.w0

    @property
    def sep(self):
        return np.array(
            [self._parameters.sep[i] for i in range(self._parameters.sep_len)]
        )

    @sep.setter
    def sep(self, value):
        self._parameters.sep = ccoffe.coffe_generate_range(0, 100, 100)
        self._parameters.sep_len = 100

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

    def compute_multipole(
        self,
        z : float,
        r : float, # should this be called r, or sep, or separations?
        l : int,
    ):
        """
        Computes the multipole of the current configuration at the point (z, r, l)
        """
        if not self._background.flag:
            self._background_init()
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

        if not self._integral.size:
            self._integrals_init()

        temp = self._parameters.output_type
        self._parameters.output_type = 1
        if not self._corrfunc.size:
            self._corrfunc_init()
        self._parameters.output_type = temp

        return np.array([
            Corrfunc(
                z=self._corrfunc.array[i].coords.z_mean,
                r=self._corrfunc.array[i].coords.separation,
                mu=self._corrfunc.array[i].coords.l,
                value=self._corrfunc.array[i].value,
            ) for i in range(self._corrfunc.size)
        ])



    # TODO implement checks to make sure that all of the structures are
    # properly initialized before running any of the below
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

