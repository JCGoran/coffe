# COFFE v3

This is the public repository for the code COFFE (COrrelation Function Full-sky Estimator), which can be used to compute the following quantities in linear perturbation theory:

* full-sky and flat-sky 2-point correlation function (2PCF) of galaxy number counts, taking into account all of the effects (density, RSD, lensing, etc.)
* full-sky and flat-sky multipoles of the 2PCF
* redshift-averaged multipoles of the 2PCF
* flat-sky Gaussian covariance matrix of the multipoles of the 2PCF
* flat-sky Gaussian covariance matrix of the redshift-averaged multipoles of the 2PCF

The relevant theoretical papers are:

* [The full-sky relativistic correlation function and power spectrum of galaxy number counts: I. Theoretical aspects, arXiv:1708.00492](https://arxiv.org/abs/1708.00492)
* [COFFE: a code for the full-sky relativistic galaxy correlation function, arXiv:1806.11090](https://arxiv.org/abs/1806.11090)
* [The flat-sky approximation to galaxy number counts - redshift space correlation function, arXiv:2011.01878](https://arxiv.org/abs/2011.01878)

## Installation

### From pip

If you are on Linux or MacOS, the latest version of COFFE can be installed using:

```sh
pip install coffe
```

Note that Windows is not officially supported.

If you wish to install the development version of COFFE, please refer to the section below.

### Development version

**NOTE**: the use of a virtual environment (such as Python's [venv](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment)) is highly recommended.

#### Prerequisites 

If you would like to install the development version, you will need to first have the following:

* a C compiler, compatible with the C99 standard
* a Python interpreter, version 3.8 or above
* [GSL](https://www.gnu.org/software/gsl/) (GNU Scientific Library) and the corresponding headers, version 2.1 or above (available as `libgsl-dev` on Debian-based, and as `gsl-devel` on RHEL/CentOS-based distros)
* [FFTW](http://www.fftw.org/download.html) and the corresponding headers, version 3 or above (available as `libfftw3-dev` on Debian-based, and as `fftw-devel` on RHEL/CentOS-based distros)
* [libconfig](http://hyperrealm.github.io/libconfig/)

Then clone this repository:

```sh
git clone https://github.com/JCGoran/coffe
```

then change directory to it:

```sh
cd coffe
```

#### Linux (CentOS/RHEL based)

Run the script:

```sh
bash scripts/install_other.sh gsl fftw libconfig
```

#### Linux (Debian/Ubuntu based)

Run the following command:

```sh
sudo apt install libgsl-dev libfftw3-dev libconfig-dev
```

and follow the instructions from the prompt.

#### MacOS (Homebrew)

You can install the necessary prerequisites using Homebrew:

```sh
brew install gsl fftw libconfig
```

**NOTE**: as a technical aside, Homebrew-installed packages (whether installed as pre-build binaries or from source) are built for the current version of your operating system. This means that they **CANNOT** be used to create a redistributable Python wheel, i.e. a wheel that works on any older version of MacOS.

#### MacOS (Conan)

As an alternative to Homebrew, one can use Conan to build the dependencies.

First install Conan using:

```sh
pip install conan
```

Then, generate a profile:

```sh
conan profile detect
```

Finally, install all of the dependencies in the `_build` directory:

```sh
conan install . --output-folder=_build --build=missing
```

Note that this may take a while as the packages are usually built from source.

**IMPORTANT NOTE**: Due to the fact that newer Apple devices have dual architectures (both `arm64` and `x86_64`), it is recommended to not mix these together, i.e. you should re-run _all_ of the above in clean `arm64` and `x86_64` environments (terminals) in separate COFFE directories to avoid any issues.

#### Installing CLASS and CUBA

COFFE also depends on the CLASS and CUBA libraries, which are not available on Homebrew or Conan, or the default Linux package repositories.
To install them, one needs to install `automake`, either via Homebrew (`brew install automake`) or via some other package manager.
They can then be built and installed by running:

```sh
bash scripts/install_other.sh class cuba
```

This will install the two packages in the directories `/opt/cuba_[ARCH]` and `/opt/class_public_[ARCH]`, where arch is either `x86_64` or `arm64` depending on your CPU architecture.

#### Installing COFFE

Now that the prerequisites are installed, you can install COFFE using:

```sh
pip install .
```

If you would additionally like to install all of the various tools for testing, generating docs, and development, you can additionally run:

```sh
pip install '.[all]'
```

## Documentation

The documentation for the latest version is available [here](https://jcgoran.github.io/coffe/).
To build the documentation, you can run `bash scripts/generate_docs.sh`.

## Bug reports and feature requests

Please use the [issue tracker](https://github.com/JCGoran/coffe/issues) to submit any bug reports and feature requests.

## License

COFFE is licensed under the GNU GPL 3.0. See the [LICENSE](LICENSE) file for more information.

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

## Testing COFFE

If you would like to test COFFE, you can do so in two ways: using either `pytest`, or `cmake` (deprecated).

### Testing with pytest

To run the tests via `pytest`, first install COFFE using the instructions above, and then run:

```sh
python -m pytest tests/
```

### Testing with `cmake` (deprecated)

If you do not want to build COFFE using `pip install`, you can instead use `cmake`, which is installable via `pip install cmake`. To do so, follow all of the above instructions, but instead of doing `pip install .`, you can instead do:

```sh
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=./install -DCOFFE_ENABLE_PYTHON=ON -DCOFFE_ENABLE_CLASS=ON -DCOFFE_ENABLE_CUBA=ON -DCOFFE_ENABLE_MATHOPTS=ON -DCOFFE_ENABLE_TESTS=ON ..
```

In case of issues with missing detection of GSL, FFTW, etc., which can happen when using Conan on MacOS, you can use:

```sh
cmake -DCMAKE_INSTALL_PREFIX=./install -DCOFFE_ENABLE_PYTHON=ON -DCOFFE_ENABLE_CLASS=ON -DCOFFE_ENABLE_CUBA=ON -DCOFFE_ENABLE_MATHOPTS=ON -DCOFFE_ENABLE_TESTS=ON -DCMAKE_TOOLCHAIN_FILE=../_build/conan_toolchain.cmake -DCOFFE_PYTHON_MINOR_VERSION=[VERSION] ..
```

where you must replace `[VERSION]` with whatever minor version of Python you are using (for instance, when using Python 3.9, replace `[VERSION]` with `9`).
Then you can build COFFE tests using:

```sh
make
```

**NOTE**: if you have Ninja installed, you can additionally pass `-G Ninja` to the above `cmake` command, and then run `ninja build` instead of `make`.

Finally, run the tests using:

```sh
ctest
```

### Building Python wheels

The building of wheels is done using the [`cibuildwheel`](https://cibuildwheel.readthedocs.io/en/stable/) utility.
To install it, run:

```sh
pip install cibuildwheel
```

#### Linux

Building of wheels on Linux requires a container engine like [Docker](https://docs.docker.com/engine/install/) or [Podman](https://podman.io/docs/installation).
Once one of those is installed, the wheels can be built using:

```sh
cibuildwheel --platform linux
```

The wheels will then be available in the `wheelhouse` subdirectory, and can then be uploaded to PyPI.

#### MacOS

The MacOS wheels require an [official Python installer](https://www.python.org/downloads/macos/); the ones from Homebrew, Conda, etc. will most likely not work.
To build the wheels, run:

```sh
cibuildwheel --platform macos
```

The wheels will then be available in the `wheelhouse` subdirectory, and can then be uploaded to PyPI.

**IMPORTANT NOTE**: if you installed GSL, FFTW, or libconfig via Brew, make sure to unlink them first using:

```sh
brew unlink gsl fftw libconfig
```

because otherwise `cibuildwheel` (or rather, `auditwheel`) may complain about mismatching OS versions.

### Releasing Python wheels

To automate the tedious task of building the wheels, they are now setup in the CI.
