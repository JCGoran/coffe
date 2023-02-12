# COFFE v3
This is the public repository for the code COFFE (COrrelation Function Full-sky Estimator), which can be used to compute the following quantities in linear perturbation theory:

* full-sky and flat-sky 2-point correlation function (2PCF) of galaxy number counts, taking into account all of the effects (density, RSD, lensing, etc.)
* full-sky and flat-sky multipoles of the 2PCF
* redshift-averaged multipoles of the 2PCF (**not available in v3**)
* flat-sky Gaussian covariance matrix of the multipoles of the 2PCF
* flat-sky Gaussian covariance matrix of the redshift-averaged multipoles of the 2PCF (**not available in v3**)

The relevant theoretical papers are:

* [The full-sky relativistic correlation function and power spectrum of galaxy number counts: I. Theoretical aspects, arXiv:1708.00492](https://arxiv.org/abs/1708.00492)
* [COFFE: a code for the full-sky relativistic galaxy correlation function, arXiv:1806.11090](https://arxiv.org/abs/1806.11090)
* [The flat-sky approximation to galaxy number counts - redshift space correlation function, arXiv:2011.01878](https://arxiv.org/abs/2011.01878)

## Installation

### From pip

If you are on Linux, the latest version of COFFE can be installed using:

```sh
pip install coffe
```

If you are on another platform, refer to the section below.

### Development version (including non-Linux machines)

If you would like to install the development version, you will need to first install the following:

* a C compiler, compatible with the C99 standard
* a Python interpreter, version 3.7 or above
* [GSL](https://www.gnu.org/software/gsl/) (GNU Scientific Library), version 2.1 or above
* [FFTW](http://www.fftw.org/download.html), version 3 or above

**NOTE**: if you are using Conda, you can install all of the above easily.
First, create a new environment:

```sh
conda create --name [NAME]
```

and activate it:

```sh
conda activate [NAME]
```

Finally, run:

```sh
conda install --channel conda-forge --file requirements.txt
```

Once you have installed the above (either natively or using Conda), you can first clone this repository:

```sh
git clone https://github.com/JCGoran/coffe
```

Then after `cd`-ing to it, run:

```sh
./install
```

to install all of the other dependencies which COFFE requires.


## Documentation

The documentation for the latest version is available [here](https://jcgoran.github.io/coffe/index.html).

## Bug reports and feature requests

Please use the [issue tracker](https://github.com/JCGoran/coffe/issues) to submit any bug reports and feature requests.
For bug reports, if you are running something other than the Docker version, please specify your platform as well as library versions.

## License

COFFE is licensed under the GNU GPL 3.0. See the `LICENSE` file for more information.

## Citations

If you use COFFE in a publication, we kindly ask that you cite the original paper describing the code, located at [arXiv:1806.11090](https://arxiv.org/abs/1806.11090).
A `bibTeX` entry is provided below for convenience.
```
@article{coffe:v1,
      author         = "Tansella, Vittorio and Jelic-Cizmek, Goran and Bonvin,
                        Camille and Durrer, Ruth",
      title          = "{COFFE: a code for the full-sky relativistic galaxy
                        correlation function}",
      year           = "2018",
      eprint         = "1806.11090",
      archivePrefix  = "arXiv",
      primaryClass   = "astro-ph.CO",
      SLACcitation   = "%%CITATION = ARXIV:1806.11090;%%"
}
```

## Development

### Testing

#### C tests (legacy)

To run the C test suite, you can use the command `make check`, which will build the binaries `test_[MODULE]`, where `[MODULE]` can currently be one of `background`, `integrals`, `corrfunc`, `multipoles`, `covariance`, and automatically run them.
Alternatively, you can build them one by one using `make test_[MODULE]`, and run them manually via `./test_[MODULE]`.
This is primarily useful when modifying the code itself, to make sure the old results of the code weren't broken by some new change (feature, bugfix, etc.).

#### Python tests

To run the Python test suite, first install the development requirements:

```sh
pip install -r pip-requirements-dev.txt
```

Then run:

```sh
pytest tests/test_coffe.py
```
