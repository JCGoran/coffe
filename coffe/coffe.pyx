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


COFFE_HUBBLE = (1./(2997.92458))


from cython.operator cimport dereference
from ctypes import CFUNCTYPE

from typing import Any, Callable, List, Tuple, Union
from dataclasses import dataclass
import os
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
        if not any([self.has_density, self.has_rsd, self.has_lensing]):
            raise ValueError(
                'No contributions specified, you need to specify at least one of \'has_density\', \'has_rsd\', \'has_lensing\''
            )


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
        _check_parameter('omega_cdm', value, (int, float), 0, 1)
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
        _check_parameter('omega_m', value, (int, float), 0, 1)
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
        _check_parameter('omega_baryon', value, (int, float), 0, 1)
        if not np.allclose(value, self.omega_baryon):
            # we set the value, rebalance the Omega budget, and free memory
            self._check_omegas(value, self.omega_cdm, self.omega_gamma)
            self._parameters.Omega0_baryon = value / self._coeff
            self._parameters.Omega0_m = self._parameters.Omega0_cdm + self._parameters.Omega0_baryon
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
        _check_parameter('omega_gamma', value, (int, float), 0, 1)
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
        _check_parameter('h', value, (int, float), 0, 1)
        if not np.allclose(value, self.h):
            # we change the big omegas first so that the product Omega * h^2 is
            # constant when we change h
            coeff = (self._parameters.h / value)**2 if self.use_little_omega else 1
            self._parameters.Omega0_m *= coeff
            self._parameters.Omega0_baryon *= coeff
            self._parameters.Omega0_gamma *= coeff
            self._parameters.Omega0_cdm = self._parameters.Omega0_m - self._parameters.Omega0_baryon
            # in order to keep things consistent in units of Mpc, we need to
            # rescale the separations when we change h
            for i in range(self._parameters.sep_len):
                self._parameters.sep[i] *= value / self._parameters.h

            # just to be on the safe side, we also need to alter the parameters
            # for the covariance
            for i in range(self._parameters.density1_len):
                self._parameters.density1[i] *= (self._parameters.h / value)**3
                self._parameters.density2[i] *= (self._parameters.h / value)**3
                self._parameters.pixelsize[i] *= value / self._parameters.h

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
        _check_parameter('w0', value, (int, float), -2, 0)
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
        _check_parameter('wa', value, (int, float), -1, 1)
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
        _check_parameter('n_s', value, (int, float), 0.5, 1.5)
        if not np.allclose(value, self.n_s):
            self._parameters.n_s = value
            self._free_except_parameters()


    @property
    def sigma8(self):
        """
        Windowed density fluctuation at r = 8 Mpc/h
        """
        return self._parameters.sigma8

    @sigma8.setter
    def sigma8(self, value):
        _check_parameter('sigma8', value, (int, float), 0, 2)
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
        _check_parameter('A_s', value, (int, float), 0, 2)
        if not np.allclose(value, self.A_s):
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
        _check_parameter('T_cmb', value, (int, float), 0, 10)
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
        _check_parameter('N_ur', value, (int, float), 0, 10)
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
        _check_parameter('m_ncdm', value, (int, float), 0, 10)
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
        _check_parameter('YHe', value, (int, float), 0, 10)
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
        _check_parameter('N_ncdm', value, int, 0, 10)
        if not np.allclose(value, self.N_ncdm):
            self._parameters.N_ncdm = value
            self._free_except_parameters()


    @property
    def sep(self):
        """
        Returns the list of separations (in Mpc) for which the
        2PCF/multipoles/covariance of multipoles should be computed.
        """
        return np.array(
            [self._parameters.sep[i] / self.h for i in range(self._parameters.sep_len)]
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
            self._parameters.sep[i] = temp[i] * self.h
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

        Note that if n > l, it returns the above multiplied by r^(n - l)
        instead (see arXiv:1806.11090, section 4.2 for further explanation).

        Parameters
        ----------
        r : float
            the separation (in Mpc) for which one wants to compute the integral

        n : int

        l : int
        """
        _check_parameter('r', r, (int, float), 0, 25000 / self.h)
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

        return ccoffe.coffe_interp_spline(&result.result, r * COFFE_HUBBLE * self.h)


    def galaxy_bias1(self, z : float):
        """
        Evaluates the galaxy bias of the first population of tracers at some
        redshift.
        """
        return evaluate_spline(&self._parameters.galaxy_bias1, z)


    def set_galaxy_bias1(self, x_sampling : List[float], y_sampling : List[float]):
        """
        Sets the value of the galaxy bias for the first population of tracers.
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
        Evaluates the galaxy bias of the second population at some redshift.
        """
        return evaluate_spline(&self._parameters.galaxy_bias2, z)


    def set_galaxy_bias2(self, x_sampling : List[float], y_sampling : List[float]):
        """
        Sets the value of the galaxy bias for the second population of tracers.
        """
        set_spline(
            &self._parameters.galaxy_bias2,
            x_sampling, y_sampling,
            self._parameters.interp_method
        )
        self._free_corrfunc()
        self._free_multipoles()


    def magnification_bias1(self, z : float):
        """
        Evaluates the magnification bias of the first population at some redshift.
        """
        return evaluate_spline(&self._parameters.magnification_bias1, z)


    def set_magnification_bias1(self, x_sampling : List[float], y_sampling : List[float]):
        """
        Sets the value of the magnification bias for the first population of
        tracers.
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
        """
        return evaluate_spline(&self._parameters.magnification_bias2, z)


    def set_magnification_bias2(self, x_sampling : List[float], y_sampling : List[float]):
        """
        Sets the value of the magnification bias of the second population of
        tracers.
        """
        set_spline(
            &self._parameters.magnification_bias2,
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
        """
        Number density of first tracers (in 1/Mpc^3) at z_mean.
        """
        return np.array(
            [self._parameters.density1[i] * self.h**3 for i in range(self._parameters.density1_len)]
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
            self._parameters.density1[i] = temp[i] / self.h**3
        self._free_corrfunc()
        self._free_multipoles()
        self._free_covariance_multipoles()


    @property
    def number_density2(self):
        """
        Number density of second tracers (in 1/Mpc^3) at z_mean.
        """
        return np.array(
            [self._parameters.density2[i] * self.h**3 for i in range(self._parameters.density2_len)]
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
            self._parameters.density2[i] = temp[i] / self.h**3
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
        """
        The pixel size of the covariance (roughly the resolution of the survey) in Mpc.
        """
        return np.array(
            [self._parameters.pixelsize[i] / self.h for i in range(self._parameters.pixelsize_len)]
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
            self._parameters.pixelsize[i] = temp[i] * self.h
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
        The half width of the redshift bin which is centered at z_mean.
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
        """
        Returns whether or not the covariance should be bin averaged (see eq. (A18) of arXiv:1509.04293).
        Note that the parameter `pixelsize` controls the bin width (in Mpc) in each redshift bin.
        If set to False, the covariance will not be bin averaged, and eq. (2.52)
        from arXiv:1806.11090 will be used.
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
        """
        Returns the maximum allowed comoving separation (in Mpc) to compute the
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
        ) / COFFE_HUBBLE / self.h


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


    def cross_spectrum(self, k : float, z1 : float, z2 : float, approximation : str = 'geometric'):
        """
        Evaluates the matter cross spectrum (in Mpc^3) at some k and z1 and z2.
        Input k must be in 1/Mpc.
        """
        _check_parameter('k', k, (int, float), self.k_min, self.k_max)
        _check_parameter('z1', z1, (int, float), 0, 15)
        _check_parameter('z2', z2, (int, float), 0, 15)

        if not self._power_spectrum_flag:
            self._power_spectrum_init()

        if not self._background.flag:
            self._background_init()

        return \
            ccoffe.coffe_interp_spline(&self._parameters.power_spectrum, k / self.h) \
           *ccoffe.coffe_interp_spline(&self._background.D1, z1) \
           *ccoffe.coffe_interp_spline(&self._background.D1, z2) / self.h**3


    def power_spectrum(self, k : float, z : float):
        """
        Evaluates the matter power spectrum at some k and z.
        Input k must be in 1/Mpc.
        """
        _check_parameter('k', k, (int, float), self.k_min, self.k_max)
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._power_spectrum_flag:
            self._power_spectrum_init()

        if not self._background.flag:
            self._background_init()

        if self._parameters.pk_type == ccoffe.COFFE_PK_LINEAR:
            return \
                ccoffe.coffe_interp_spline(&self._parameters.power_spectrum, k / self.h) \
               *ccoffe.coffe_interp_spline(&self._background.D1, z)**2 / self.h**3
        return ccoffe.coffe_interp_spline2d(&self._parameters.power_spectrum2d, z, k / self.h) / self.h**3


    @property
    def k_min(self):
        """
        The minimum wavenumber (in 1/Mpc) for which the power spectrum should be computed.
        """
        return self._parameters.k_min * self.h

    @k_min.setter
    def k_min(self, value):
        _check_parameter('k_min', value, (int, float), 1e-7 * self.h, 1e-3 * self.h)
        self._parameters.k_min = value / self.h
        self._parameters.k_min_norm = value / COFFE_HUBBLE / self.h


    @property
    def k_max(self):
        """
        The maximum wavenumber (in 1/Mpc) for which the power spectrum should be computed.
        """
        return self._parameters.k_max * self.h

    @k_max.setter
    def k_max(self, value):
        _check_parameter('k_max', value, (int, float), 0.1 * self.h, 1e3 * self.h)
        self._parameters.k_max = value / self.h
        self._parameters.k_max_norm = value / COFFE_HUBBLE / self.h


    def set_power_spectrum_linear(self, k : List[float], pk : List[float], z : float = 0):
        """
        Sets the linear matter power spectrum, optionally at some redshift (by default it's
        assumed the input is at z = 0).
        The input k must be in units 1/Mpc, and P(k) in units Mpc^3.

        Parameters
        ----------
        k : List[float]
            the list of wavenumbers

        pk : List[float]
            the list of values of the linear power spectrum

        z : float, default = 0
            the redshift at which the power spectrum is evaluated
        """
        # value checking
        if not np.allclose(k, np.sort(k)):
            raise ValueError('The input k array must be sorted')

        if len(k) != len(pk):
            raise ValueError('The input arrays have mismatching lengths')

        # TODO error checking for z and rescaling if z > 0

        self._free_power_spectrum()
        self._parameters.k_min = k[0] / self.h
        self._parameters.k_max = k[-1] / self.h
        self._parameters.k_min_norm = k[0] / COFFE_HUBBLE / self.h
        self._parameters.k_max_norm = k[-1] / COFFE_HUBBLE / self.h
        size = len(k)

        cdef double *x = <double *>ccoffe.coffe_malloc(sizeof(double) * size)
        cdef double *y = <double *>ccoffe.coffe_malloc(sizeof(double) * size)
        cdef double *x_norm = <double *>ccoffe.coffe_malloc(sizeof(double) * size)
        cdef double *y_norm = <double *>ccoffe.coffe_malloc(sizeof(double) * size)

        for (i, ki), pki in zip(enumerate(k), pk):
            x[i] = ki / self.h
            y[i] = pki * self.h**3
            x_norm[i] = ki / COFFE_HUBBLE / self.h
            y_norm[i] = pki * COFFE_HUBBLE**3 * self.h**3

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
        """
        Evaluates the comoving distance at some redshift (in Mpc).
        """
        # TODO figure out how to dereference an operator
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._background.flag:
            self._background_init()

        return ccoffe.coffe_interp_spline(&self._background.comoving_distance, z) / COFFE_HUBBLE / self.h


    def hubble_rate(self, z : float):
        """
        Evaluates the Hubble rate (H(z)) at some redshift (in 1/Mpc).
        """
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._background.flag:
            self._background_init()

        return ccoffe.coffe_interp_spline(&self._background.Hz, z) * COFFE_HUBBLE * self.h


    def scale_factor(self, z : float):
        """
        Returns the scale factor a evaluated at some redshift.
        """
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._background.flag:
            self._background_init()

        return ccoffe.coffe_interp_spline(&self._background.a, z)


    def hubble_rate_conformal(self, z : float):
        """
        Evaluates the conformal Hubble rate (𝓗(z)) at some redshift (in 1/Mpc).
        """
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._background.flag:
            self._background_init()

        return ccoffe.coffe_interp_spline(&self._background.conformal_Hz, z) * COFFE_HUBBLE * self.h


    def hubble_rate_conformal_derivative(self, z : float):
        """
        Evaluates the first derivative of the conformal Hubble rate w.r.t.
        conformal time (d𝓗(z)/dτ) at some redshift (in 1/Mpc^2).
        """
        _check_parameter('z', z, (int, float), 0, 15)

        if not self._background.flag:
            self._background_init()

        return ccoffe.coffe_interp_spline(&self._background.conformal_Hz_prime, z) * COFFE_HUBBLE**2 * self.h**2


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
        recompute : bool = False,
    ):
        """
        Computes the correlation function of the current configuration at the
        point (z, r, mu)

        Parameters
        ----------
        z : float
            the mean redshift

        r : float
            the comoving separation (in Mpc) between the two points in the sky

        mu : float
            the angle mu (see `help(coffe.mu)` for definition)
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
        _check_parameter('r', r, (int, float), 0, 1500 / self.h)
        _check_parameter('mu', mu, (int, float), -1, 1)

        cdef ccoffe.gsl_error_handler_t *default_handler = \
            ccoffe.gsl_set_error_handler_off()

        result = ccoffe.coffe_integrate(
            &self._parameters,
            &self._background,
            &self._integral,
            z, r * COFFE_HUBBLE * self.h, mu, 0,
            ccoffe.NONINTEGRATED, ccoffe.CORRFUNC
        ) + \
        ccoffe.coffe_integrate(
            &self._parameters,
            &self._background,
            &self._integral,
            z, r * COFFE_HUBBLE * self.h, mu, 0,
            ccoffe.SINGLE_INTEGRATED, ccoffe.CORRFUNC
        ) + \
        ccoffe.coffe_integrate(
            &self._parameters,
            &self._background,
            &self._integral,
            z, r * COFFE_HUBBLE * self.h, mu, 0,
            ccoffe.DOUBLE_INTEGRATED, ccoffe.CORRFUNC
        )

        ccoffe.gsl_set_error_handler(default_handler)

        return result


    def compute_multipole(
        self, *,
        z : float, r : float, l : int,
        recompute : bool = False,
    ):
        """
        Computes the multipole of the current configuration at the point (z, r, l)

        Parameters
        ----------
        z : float
            the mean redshift

        r : float
            the comoving separation (in Mpc) between the 2 points in the sky

        l : int
            the multipole moment
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
        _check_parameter('r', r, (int, float), 0, 1500 / self.h)
        _check_parameter('l', l, (int,), 0, 10)

        cdef ccoffe.gsl_error_handler_t *default_handler = \
            ccoffe.gsl_set_error_handler_off()

        result = ccoffe.coffe_integrate(
            &self._parameters,
            &self._background,
            &self._integral,
            z, r * COFFE_HUBBLE * self.h, 0, l,
            ccoffe.NONINTEGRATED, ccoffe.MULTIPOLES
        ) + \
        ccoffe.coffe_integrate(
            &self._parameters,
            &self._background,
            &self._integral,
            z, r * COFFE_HUBBLE * self.h, 0, l,
            ccoffe.SINGLE_INTEGRATED, ccoffe.MULTIPOLES
        ) + \
        ccoffe.coffe_integrate(
            &self._parameters,
            &self._background,
            &self._integral,
            z, r * COFFE_HUBBLE * self.h, 0, l,
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
        an array of instances of `Multipoles`.
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
                r=self._multipoles.array[i].coords.separation / COFFE_HUBBLE / self.h,
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
        an array of instances of `Corrfunc`.

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
                r=self._corrfunc.array[i].coords.separation / COFFE_HUBBLE / self.h,
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
        an array of instances of `Covariance`.
        """
        self._check_coords_multipoles()
        self._check_contributions()
        recompute = bool(recompute)

        if not all(
            len(self.z_mean) == len(_) \
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
                r1=self._covariance_multipoles.array[i].coords.separation1 / self.h,
                r2=self._covariance_multipoles.array[i].coords.separation2 / self.h,
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
