import os
import numpy as np
import pandas as pd
import pytest

import coffe

DATA_DIR = 'tests/benchmarks/'

class TestCoffe:
    def test_bias(self):
        """
        Tests for setting and getting the biases.
        """
        cosmo = coffe.Coffe()

        with pytest.raises(ValueError):
            x = np.linspace(-1, 10, 100)
            y = np.zeros(len(x))
            cosmo.set_galaxy_bias1(x, y)

        with pytest.raises(ValueError):
            x = np.linspace(0, 100, 100)
            y = np.zeros(len(x))
            cosmo.set_galaxy_bias1(x, y)

        with pytest.raises(ValueError):
            x = np.linspace(0, 1, 100)
            y = np.zeros(len(x))
            cosmo.set_galaxy_bias1(x, y)
            cosmo.galaxy_bias1(2)


    def test_parameters(self):
        """
        Tests setting of various cosmological parameters.
        """
        cosmo = coffe.Coffe()

        with pytest.raises(ValueError):
            cosmo.omega_cdm = 1.1

        with pytest.raises(ValueError):
            cosmo.h = -1



    def test_background(self):
        cosmo = coffe.Coffe()
        # we can't init the background automatically, so this is cheating a bit
        cosmo.comoving_distance(1)

        data = np.loadtxt(os.path.join(DATA_DIR, 'benchmark_background.dat'))

        d = {name : data[:, index] \
            for index, name in enumerate([
                'z', 'scale_factor', 'hubble_rate', 'hubble_rate_conformal', \
                'hubble_rate_conformal_prime', 'growth_factor', 'growth_rate', \
                'conformal_distance'
            ])
        }

        for key in d:
            if hasattr(cosmo, key):
                # I forgot to normalize the growth factor in the test
                if key == 'growth_factor':
                    for index, zi in enumerate(d['z']):
                        assert np.isclose(d[key][index] / d[key][0], getattr(cosmo, key)(zi))
                else:
                    for index, zi in enumerate(d['z']):
                        assert np.isclose(d[key][index], getattr(cosmo, key)(zi))


    def test_integrals(self):
        cosmo = coffe.Coffe(
            # in Mpc/h
            sep=[10, 20, 40, 100, 150],
            mu=[0.0, 0.2, 0.5, 0.8, 0.95],
        )

        k, pk = np.transpose(np.loadtxt('PkL_CLASS.dat'))
        cosmo.set_power_spectrum_linear(k, pk)

        # mapping of indices to (n, l) pairs
        mapping = {
            0 : {'l' : 0, 'n' : 0},
            1 : {'l' : 2, 'n' : 0},
            2 : {'l' : 4, 'n' : 0},
            3 : {'l' : 1, 'n' : 1},
            4 : {'l' : 3, 'n' : 1},
            5 : {'l' : 0, 'n' : 2},
            6 : {'n' : 2, 'l' : 2},
            7 : {'n' : 3, 'l' : 1},
        }

        for index in mapping:
            data = np.loadtxt(os.path.join(DATA_DIR, f'benchmark_integral{index}.dat'))
            xarr, yarr = np.transpose(data)
            for x, y in zip(xarr, yarr):
                if x / coffe._COFFE_HUBBLE > 20000:
                    assert np.isclose(
                        y,
                        cosmo.integral(
                            x / coffe._COFFE_HUBBLE,
                            l=mapping[index]['l'],
                            n=mapping[index]['n'],
                        )
                    )


    def test_corrfunc(self):
        cosmo = coffe.Coffe(
            # in Mpc/h
            sep=[10, 20, 40, 100, 150],
            mu=[0.0, 0.2, 0.5, 0.8, 0.95],
        )

        k, pk = np.transpose(np.loadtxt('PkL_CLASS.dat'))

        contributions = {'den' : 'density', 'rsd' : 'rsd', 'len' : 'lensing'}
        for prefix in contributions:
            cosmo.reset_contributions()
            setattr(cosmo, f'has_{contributions[prefix]}', True)
            cosmo.set_power_spectrum_linear(k, pk)
            result = cosmo.compute_corrfunc_bulk()
            df = pd.DataFrame([_.to_dict() for _ in result])
            for index, mu in enumerate(cosmo.mu):
                data = np.loadtxt(
                    os.path.join(DATA_DIR, f'benchmark_{prefix}_corrfunc{index}.dat')
                )
                x, y = np.transpose(data)
                assert np.allclose(df.loc[df.mu == mu].r.values, x)
                assert np.allclose(df.loc[df.mu == mu].value.values, y, rtol=5e-4)


    def test_multipoles(self):
        cosmo = coffe.Coffe(
            # in Mpc/h
            sep=[10, 20, 40, 100, 150],
            l=[0, 2, 4],
        )

        k, pk = np.transpose(np.loadtxt('PkL_CLASS.dat'))
        cosmo.set_power_spectrum_linear(k, pk)

        contributions = {'den' : 'density', 'rsd' : 'rsd', 'len' : 'lensing'}
        for prefix in contributions:
            cosmo.reset_contributions()
            setattr(cosmo, f'has_{contributions[prefix]}', True)
            cosmo.set_power_spectrum_linear(k, pk)
            result = cosmo.compute_multipoles_bulk()
            df = pd.DataFrame([_.to_dict() for _ in result])
            for mp in cosmo.l:
                data = np.loadtxt(
                    os.path.join(DATA_DIR, f'benchmark_{prefix}_multipoles{mp}.dat')
                )
                x, y = np.transpose(data)
                assert np.allclose(df.loc[df.l == mp].r.values, x)
                assert np.allclose(df.loc[df.l == mp].value.values, y, rtol=5e-4)


    def test_covariance_multipoles(self):
        cosmo = coffe.Coffe()
        cosmo.set_parameters(
            has_density=True,
            has_rsd=True,
            number_density=[1e-3],
            z_mean=[1.0],
            deltaz=[0.1],
            fsky=[0.2],
            pixelsize=[50],
            sep=np.arange(50, 350, 50),
            l=[0, 2, 4],
        )

        k, pk = np.transpose(np.loadtxt('PkL_CLASS.dat'))
        cosmo.set_power_spectrum_linear(k, pk)

        result = cosmo.compute_covariance_bulk()
        df = pd.DataFrame([_.to_dict() for _ in result])
        for mp1 in cosmo.l:
            for mp2 in cosmo.l:
                data = np.loadtxt(
                    os.path.join(DATA_DIR, f'benchmark_multipoles_covariance_{mp1}{mp2}.dat')
                )
                x, y, z = np.transpose(data)
                assert np.allclose(df.loc[(df.l1 == mp1) & (df.l2 == mp2)].r1.values, x)
                assert np.allclose(df.loc[(df.l1 == mp1) & (df.l2 == mp2)].r2.values, y)
                assert np.allclose(df.loc[(df.l1 == mp1) & (df.l2 == mp2)].value.values, z, rtol=5e-4)

        assert np.allclose(
            cosmo.covariance_matrix(),
            np.transpose(cosmo.covariance_matrix())
        )

        assert np.allclose(
            cosmo.covariance_matrix_inverse(),
            np.linalg.inv(cosmo.covariance_matrix())
        )
