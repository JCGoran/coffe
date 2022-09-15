import os
import numpy as np
import pandas as pd
import pytest

import coffe
from coffe_utils import covariance_matrix

DATA_DIR = "tests/benchmarks/"


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

        # should raise since omega_baryon = 0.05 by default
        with pytest.raises(ValueError):
            cosmo.omega_cdm = 0.99

        with pytest.raises(ValueError):
            cosmo.omega_m = 0.99999999

        with pytest.raises(ValueError):
            cosmo.z_mean = (-1,)

        with pytest.raises(TypeError):
            cosmo.z_mean = 0.1

        with pytest.raises(ValueError):
            cosmo.z_mean = []
            cosmo.compute_multipoles_bulk()

        with pytest.raises(ValueError):
            cosmo.z_mean = []
            cosmo.compute_corrfunc_bulk()

        with pytest.raises(ValueError):
            cosmo.z_mean = []
            cosmo.compute_covariance_bulk()

        with pytest.raises(ValueError):
            cosmo.reset_contributions()
            cosmo.compute_corrfunc_bulk()

        # still has no contributions
        with pytest.raises(ValueError):
            cosmo.compute_corrfunc(z=0.1, r=100, mu=0.5)

        with pytest.raises(ValueError):
            cosmo.compute_multipole(z=0.1, r=100, l=2)

        with pytest.raises(ValueError):
            cosmo.integral(r=100, n=4, l=4)

    def test_background(self):
        cosmo = coffe.Coffe()
        # we can't init the background automatically, so this is cheating a bit
        cosmo.comoving_distance(1)

        data = np.loadtxt(os.path.join(DATA_DIR, "benchmark_background.dat"))

        d = {
            name: data[:, index]
            for index, name in enumerate(
                [
                    "z",
                    "scale_factor",
                    "hubble_rate",
                    "hubble_rate_conformal",
                    "hubble_rate_conformal_prime",
                    "growth_factor",
                    "growth_rate",
                    "conformal_distance",
                ]
            )
        }

        for key in d:
            if hasattr(cosmo, key):
                # I forgot to normalize the growth factor in the test
                if key == "growth_factor":
                    for index, zi in enumerate(d["z"]):
                        assert np.isclose(
                            d[key][index] / d[key][0], getattr(cosmo, key)(zi)
                        )
                else:
                    for index, zi in enumerate(d["z"]):
                        assert np.isclose(d[key][index], getattr(cosmo, key)(zi))

    def test_integrals(self):
        cosmo = coffe.Coffe(
            # in Mpc/h
            sep=[10, 20, 40, 100, 150],
            mu=[0.0, 0.2, 0.5, 0.8, 0.95],
        )

        k, pk = np.transpose(np.loadtxt("PkL_CLASS.dat"))
        cosmo.set_power_spectrum_linear(k, pk)

        # mapping of indices to (n, l) pairs
        mapping = {
            0: {"l": 0, "n": 0},
            1: {"l": 2, "n": 0},
            2: {"l": 4, "n": 0},
            3: {"l": 1, "n": 1},
            4: {"l": 3, "n": 1},
            5: {"l": 0, "n": 2},
            6: {"n": 2, "l": 2},
            7: {"n": 3, "l": 1},
        }

        for index in mapping:
            data = np.loadtxt(os.path.join(DATA_DIR, f"benchmark_integral{index}.dat"))
            xarr, yarr = np.transpose(data)
            for x, y in zip(xarr, yarr):
                if x / coffe.COFFE_HUBBLE > 1 and x / coffe.COFFE_HUBBLE < 20000:
                    assert np.isclose(
                        y,
                        cosmo.integral(
                            r=x / coffe.COFFE_HUBBLE,
                            l=mapping[index]["l"],
                            n=mapping[index]["n"],
                        ),
                    )

    def test_corrfunc(self):
        cosmo = coffe.Coffe(
            # in Mpc/h
            sep=[10, 20, 40, 100, 150],
            mu=[0.0, 0.2, 0.5, 0.8, 0.95],
        )

        k, pk = np.transpose(np.loadtxt("PkL_CLASS.dat"))

        contributions = {
            "den": "density",
            "rsd": "rsd",
            "len": "lensing",
            "d1": "d1",
            "d2" : "d2",
            "g1": "g1",
            "g2": "g2",
            "g3": "g3",
        }
        for prefix in contributions:
            cosmo.reset_contributions()
            setattr(cosmo, f"has_{contributions[prefix]}", True)
            cosmo.set_power_spectrum_linear(k, pk)
            result = cosmo.compute_corrfunc_bulk()
            df = pd.DataFrame([_.to_dict() for _ in result])
            for index, mu in enumerate(cosmo.mu):
                data = np.loadtxt(
                    os.path.join(DATA_DIR, f"benchmark_{prefix}_corrfunc{index}.dat")
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

        k, pk = np.transpose(np.loadtxt("PkL_CLASS.dat"))
        cosmo.set_power_spectrum_linear(k, pk)

        contributions = {"den": "density", "rsd": "rsd", "len": "lensing"}
        for prefix in contributions:
            cosmo.reset_contributions()
            setattr(cosmo, f"has_{contributions[prefix]}", True)
            cosmo.set_power_spectrum_linear(k, pk)
            result = cosmo.compute_multipoles_bulk()
            df = pd.DataFrame([_.to_dict() for _ in result])
            for mp in cosmo.l:
                data = np.loadtxt(
                    os.path.join(DATA_DIR, f"benchmark_{prefix}_multipoles{mp}.dat")
                )
                x, y = np.transpose(data)
                assert np.allclose(df.loc[df.l == mp].r.values, x)
                assert np.allclose(df.loc[df.l == mp].value.values, y, rtol=5e-4)

    def test_covariance_multipoles(self):
        cosmo = coffe.Coffe()
        cosmo.set_parameters(
            has_density=True,
            has_rsd=True,
            number_density1=[1e-3],
            number_density2=[1e-3],
            z_mean=[1.0],
            deltaz=[0.1],
            fsky=[0.2],
            pixelsize=[50],
            sep=np.arange(50, 350, 50),
            l=[0, 2, 4],
        )

        k, pk = np.transpose(np.loadtxt("PkL_CLASS.dat"))
        cosmo.set_power_spectrum_linear(k, pk)

        result = cosmo.compute_covariance_bulk()
        df = pd.DataFrame([_.to_dict() for _ in result])
        for mp1 in cosmo.l:
            for mp2 in cosmo.l:
                data = np.loadtxt(
                    os.path.join(
                        DATA_DIR, f"benchmark_multipoles_covariance_{mp1}{mp2}.dat"
                    )
                )
                x, y, z = np.transpose(data)
                assert np.allclose(df.loc[(df.l1 == mp1) & (df.l2 == mp2)].r1.values, x)
                assert np.allclose(df.loc[(df.l1 == mp1) & (df.l2 == mp2)].r2.values, y)
                assert np.allclose(
                    df.loc[(df.l1 == mp1) & (df.l2 == mp2)].value.values, z, rtol=5e-4
                )

        assert covariance_matrix(result).ndim == 2
        assert np.shape(covariance_matrix(result)) == (18, 18)

        assert np.allclose(
            covariance_matrix(result), np.transpose(covariance_matrix(result))
        )

        # we don't need to test the trace or the determinant since that follows
        # automatically if the matrix has all positive eigenvalues
        assert np.all(np.linalg.eigvalsh(covariance_matrix(result)) > 0)

        assert covariance_matrix(result, rstep=111).size == 0

    def test_error_handler(self):
        """
        Checks that the GSL error handler doesn't abort the computation
        """
        # initialize COFFE with some cosmology
        cosmo = coffe.Coffe()
        cosmo.set_parameters(
            sep=np.linspace(10, 300, 100),
            has_rsd=False,
            omega_m=0.31,
        )
        cosmo.compute_multipole(l=4, r=10, z=1.5)

    def test_multiple_populations(self):
        """
        Tests for 2 populations of galaxies
        """
        # no covariance contributions
        cosmo = coffe.Coffe(
            has_density=True,
            has_rsd=True,
            number_density1=[1e-5],
            number_density2=[1e-5],
            z_mean=[1.0],
            deltaz=[0.1],
            fsky=[0.2],
            pixelsize=[50],
            sep=np.arange(50, 350, 50),
            l=[0, 2, 4],
            covariance_cosmic=False,
            covariance_mixed=False,
            covariance_poisson=False,
        )

        k, pk = np.transpose(np.loadtxt("PkL_CLASS.dat"))
        cosmo.set_power_spectrum_linear(k, pk)

        cov = covariance_matrix(cosmo.compute_covariance_bulk())

        assert np.allclose(cov, 0)

        # only Poisson contribution
        cosmo.covariance_poisson = True

        k, pk = np.transpose(np.loadtxt("PkL_CLASS.dat"))
        cosmo.set_power_spectrum_linear(k, pk)

        cov = covariance_matrix(cosmo.compute_covariance_bulk())

        # the diagonal should not be zero
        assert not np.allclose(np.diag(cov), 0)

        # for 4 different poulations, covariance with mixed and Poisson
        # contributions only should be zero
        cosmo.set_parameters(
            covariance_populations=[1, 2, 3, 4],
            covariance_poisson=True,
            covariance_mixed=True,
            covariance_cosmic=False,
        )

        cov = covariance_matrix(cosmo.compute_covariance_bulk())

        assert np.allclose(cov, 0)


class TestRepresentation:
    """
    Tests for Corrfunc, Multipoles, and Covariance classes.
    """

    def test_representation(self):
        with pytest.raises(ImportError):
            from coffe import Representation

        coffe.Multipoles(l=0, r=10, z=1.0, value=1e-3)
