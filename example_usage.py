# Example usage of the Python wrapper
import numpy as np
import coffe



# setting the parameters
cosmology = coffe.Coffe(
    **{
        'z_mean' : [1.5],
        'has_density' : True,
        'has_rsd' : True,
        'mu' : [0.0, 0.1, 0.5, 0.7, 0.9, 0.95, 0.98, 0.99],
        'omega_baryon' : 0.04,
        'omega_cdm' : 0.26,
        'sep' : np.loadtxt('separations.dat'),
    },
)

# computes the 2PCF for the above parameters
result = cosmology.compute_corrfunc_bulk()
print(result)
