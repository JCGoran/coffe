# Example usage of the Python wrapper
# note that the COFFE binary must be already built
# using `make`
# to run this script itself, use `python example_usage.py`
import sys

# this is needed to correctly import COFFE (from the current directory)
sys.path.append('python/')

import coffe

# running with some parameters
coffe.run_coffe(
    parameters={
        'z_mean' : 1.5,
        'correlation_contributions' : ['den', 'rsd'],
        'mu' : [0.0, 0.1, 0.5, 0.7, 0.9, 0.95, 0.98, 0.99],
        'omega_baryon' : 0.04,
        'omega_cdm' : 0.26,
        'output_path' : 'python_wrapper_test/',
        'output_prefix' : '$TIME',
        'output_type' : 1,
        'input_separations' : 'separations.dat',
    },
    cores=1,
    path='./coffe',
)
