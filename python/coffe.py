"""
Simple Python wrapper for COFFE
"""

import os

import tempfile
import subprocess as sp
import libconf as lc

# the default parameters to be used by COFFE
def get_default_parameters():
    """
    Returns the default parameters as used by COFFE
    in the form of a Python dictionary.
    """
    return {
        'input_separations' : 'separations.dat',
        'input_power_spectrum' : 'PkL_CLASS.dat',
        'omega_cdm' : 0.25,
        'omega_baryon' : 0.05,
        'omega_gamma' : 9e-5,
        'w0' : -1.0,
        'wa' : 0.0,
        'galaxy_bias1' : 1.,
        'read_galaxy_bias1' : 0,
        'input_galaxy_bias1' : '',
        'galaxy_bias2' : 1.,
        'galaxy_bias_analytic' : 0,
        'read_galaxy_bias2' : 0,
        'input_galaxy_bias2' : '',
        'magnification_bias1' : 0.,
        'read_magnification_bias1' : 0,
        'input_magnification_bias1' : '',
        'magnification_bias2' : 0.,
        'read_magnification_bias2' : 0,
        'input_magnification_bias2' : '',
        'evolution_bias1' : 0.,
        'read_evolution_bias1' : 0,
        'input_evolution_bias1' : '',
        'evolution_bias2' : 0.,
        'read_evolution_bias2' : 0,
        'input_evolution_bias2' : '',
        'covariance_z_mean' : [0.5, 1.0, 1.5],
        'covariance_deltaz' : [0.1, 0.1, 0.2],
        'covariance_density' : [0.1, 0.2, 0.01],
        'covariance_fsky' : [0.3, 0.5, 0.2],
        'covariance_pixelsize' : 10.0,
        'covariance_zmin' : [2.0, 2.2, 2.3],
        'covariance_zmax' : [2.5, 2.8, 2.5],
        'output_path' : 'results/',
        'output_prefix' : '$TIME',
        'correlation_contributions' : ['den'],
        'output_type' : 2,
        'z_mean' : 1.0,
        'deltaz' : 0.1,
        'z_min' : 0.7,
        'z_max' : 1.3,
        'mu' : [0.7],
        'multipoles' : [0, 2, 4],
        'output_background' : [
            'z', 'a', 'H',
            'conformal_H', 'conformal_H_prime',
            'D1', 'f', 'comoving_distance',
        ],
        'flatsky_local' : 0,
        'flatsky_density_lensing' : 0,
        'flatsky_lensing_lensing' : 0,
        'background_sampling' : 10000,
        'bessel_sampling' : 10000,
        'theta_sampling' : 3000,
        'integration_method' : 2,
        'integration_sampling' : 750000,
        'k_min' : 1e-5,
        'k_max' : 300.,
        'interpolation' : 5,
        'have_class' : 1,
        'pk_type' : 0,
        'zeldovich_approximation' : 0,
        'h' : 0.67,
        'k_pivot' : 0.05,
        'ln_10_pow_10_A_s' : 3.06,
        'n_s' : 0.96,
        'have_window' : 0,
        'window_size' : 5,
        'covariance_integration_method' : 1,
        'covariance_integration_bins' : 8000,
        'covariance_interpolation_method' : 2,
        'covariance_minimum_separation' : 20,
        'verbose' : 1,
        'only_cross_correlations' : 0,
    }

def _check_value(value):
    if value:
        try:
            os.stat(value)
        except:
            raise ValueError(
                f'File \'{value}\' does not exist!'
            )

# actually running COFFE
def run_coffe(
        parameters={},
        cores=8,
        path='./coffe',
):
    """
    Semi-useful wrapper for running COFFE from Python.
    Just set the path to wherever the COFFE binary is,
    and give a bunch of parameters as input in the form
    of a Python dictionary.
    """

    # check COFFE is where the user says it is
    try:
        os.stat(path)
    except:
        raise ValueError(
            f'The COFFE binary \'{path}\' does not exist,' \
            'make sure the `path` variable is set correctly!'
        )

    # overwrite defaults if input is something
    default_parameters = {
        **get_default_parameters(),
        **parameters,
    } if isinstance(parameters, dict) else get_default_parameters()

    for _name in ['matter', 'magnification', 'evolution']:
        for _index in [1, 2]:
            if default_parameters.get(f'read_{_name}_bias{_index}'):
                _check_value(default_parameters.get(f'input_{_name}_bias{_index}'))

    settingsfile = tempfile.NamedTemporaryFile(
        mode='w',
        delete=False,
    )

    # put settings into a temporary settingsfile
    settingsfile.write(
        lc.dumps(default_parameters)
    )
    settingsfile.close()

    # run COFFE
    sp.run(
        [
            f'{path}',
            '-s', f'{settingsfile.name}',
            '-n', f'{cores}',
        ],
        check=True,
    )

    # garbage collection
    os.unlink(settingsfile.name)
