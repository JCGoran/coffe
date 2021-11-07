# distutils: sources = src/errors.c src/common.c src/parser.c src/background.c
# distutils: include_dirs = src/ ./
# distutils: libraries = m gsl gslcblas config

cimport ccoffe

import numpy as np

cdef class Coffe:
    cdef ccoffe.coffe_parameters_t _parameters
    cdef ccoffe.coffe_background_t _background
    cdef ccoffe.coffe_corrfunc_array_t _corrfunc

    def __cinit__(self):
        ccoffe.coffe_parse_default_parameters(&self._parameters)
        #self._background.flag = 0
        self._corrfunc.size = 0
        self._corrfunc.array = NULL

    def __dealloc__(self):
        ccoffe.coffe_parameters_free(&self._parameters)

    def _balance_content(self):
        """
        Balances the energy content of dark energy so it all adds up to 1
        """
        self._parameters.Omega0_de = 1 - self._parameters.Omega0_cdm - self._parameters.Omega0_baryon - self._parameters.Omega0_gamma

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
