#!/usr/bin/env python3

# Example usage of the Python wrapper

import numpy as np

import coffe

if __name__ == "__main__":
    # setting the parameters
    cosmology = coffe.Coffe(
        **{
            "z_mean": [1.5],
            "deltaz": [0.1],
            "has_density": True,
            "has_rsd": True,
            "mu": [0.0, 0.1, 0.5, 0.7, 0.9, 0.95, 0.98, 0.99],
            "omega_baryon": 0.04,
            "omega_m": 0.31,
            "sep": np.arange(5, 100, 5),
        },
    )

    # computes the 2PCF for the above parameters
    result = cosmology.compute_corrfunc_bulk()
    print(result)
