import os
import numpy as np
from scipy.interpolate import interp1d
import pytest

import coffe

DATA_DIR = 'tests/benchmarks/'

class TestCoffe:
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
