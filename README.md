# COFFE
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

## Installation and running
There are a couple of ways to install and run the code:

### 1. From the tarball
Go [here](https://github.com/JCGoran/coffe/releases) and download the latest release that doesn't have "Docker" in the title.
Extract the contents of the tarball somewhere.
Make sure you have a C compiler (compatible with the C99 standard) and the following libraries installed:
* [FFTW](http://www.fftw.org/download.html)
* [libconfig](https://hyperrealm.github.io/libconfig/)
* GSL
* [CUBA](http://www.feynarts.de/cuba/) (optional)
* [CLASS](https://github.com/lesgourg/class_public) (optional)

Once you have installed the dependencies, you can run 
```
./configure && make
```
from the directory in which COFFE was extracted.
If you wish to enable the CUBA library, you can specify the flag `--enable-cuba` for the `configure` script.
You can also enable CLASS as a library using `--enable-class` and double exponential quadrature using `--enable-doubleexp`.

The CLASS and CUBA libraries as well as their corresponding header files need to be installed in the search path of the compiler and the linker, otherwise you **must** manually set them, by running `configure` with additional options which specify the paths:

```
./configure \
    CPPFLAGS="$CPPFLAGS -I/path/to/class_headers/ -I/path/to/cuba_headers/" \
    LDLAGS="$LDFLAGS -L/path/to/class_library/ -L/path/to/cuba_library/" \
    LIBS="$LIBS -lclass -lcuba"
```

You can get the default search path for the headers using the command (tested with `gcc`):

```
cc -v -E -xc /dev/null -o /dev/null
```

and the corresponding one for the libraries (using `ld`):

```
ld --verbose | grep SEARCH_DIR | tr -s ' ;' \\012
```

To run the code:
```
coffe -s [SETTINGSFILE] -n [NUMTHREADS]
```
where `[SETTINGSFILE]` is the name of the settings file, and `[NUMTHREADS]` is the number of processes you wish to spawn (loosely speaking, the number of cores to use for the computation).

The `settings.cfg` file contains explanations about the possible input and output.
For more details, please consult the manual located in the `manual` subdirectory.

### 2. From source
After making sure you have the adequate versions of `autotools` installed on your machine, clone this repository:
```
git clone https://github.com/JCGoran/coffe
```
and run
```
autoreconf -i
```
From there, you may proceed as under point 1.

### 3. From DockerHub (compatible with Docker, Singularity and Shifter)
COFFE has a [Docker](https://docs.docker.com/install/) version available, which comes packaged with the CLASS and CUBA libraries.
To pull the Docker image from DockerHub:
```
docker pull jcgoran/coffe:[TAG]
```
and then run a container using
```
docker run -ti jcgoran/coffe:[TAG]
```
where `[TAG]` is the tag of the release (`latest` by default).
The `coffe` binary can be run from any directory _inside_ the container.
Copying the input/output can be done using `docker cp` as described in the manual.
To run the code without having to use `docker cp`, you may run COFFE as follows:
```
docker run --rm -ti -v "$(pwd):/data" jcgoran/coffe:[TAG] coffe -s [SETTINGSFILE] -n [NUMTHREADS]
```
This will automatically mount your current directory onto the container, and COFFE will automatically write the output to your storage device, just like the version built from source.
Note that specifying a settings file in a parent directory (i.e. using `../` and variations of it) or using absolute paths (i.e. `/path/to/dir` instead of `path/to/dir`) will not work when using the above command.
As the above command is a bit lengthy, you could put the following alias in your `.bashrc` file:
```
alias coffe='docker run --rm -ti -v "$(pwd):/data" jcgoran/coffe:[TAG] coffe'
```
This allows you to run COFFE from any directory by simply calling:
```
coffe -s [SETTINGSFILE] -n [NUMTHREADS]
```
**NOTE**: when running the Docker version and using `output=$TIME` in the settings file, the system time of the Docker image is not necessarily synchronized with your local time.
To circumvent this, you may run the above command with an additional mount flag `-v /etc/localtime:/etc/localtime:ro`, which should synchronize the system time of the Docker image with your local time.

**NOTE 2**: due to Docker's permission model, the default permissions of the output files may be wrong when using the `-v [MOUNT]` flag, i.e. they may belong to `root`, so you may not be able to access them; to fix this, you can either:

1. run `sudo chmod -R $USER:$USER [DIRECTORY]` on the output directory every you run the Docker image
2. add `--user "$(id -u):$(id -g)"` to the above Docker command

### Docker alternatives
Alternatively, one can run the Docker image using [Singularity](https://github.com/sylabs/singularity).
To pull the image from DockerHub using Singularity:
```
singularity pull docker://jcgoran/coffe
```
and then run the code directly using
```
singularity exec [FILE] coffe -s [SETTINGSFILE] -n [NUMTHREADS]
```
where `[FILE]` is the name of the file Singularity has downloaded.
This, unlike Docker, avoids the hassle of copying to/from the container and directly outputs to your storage device.

Alternatively, one can use [Shifter](https://github.com/NERSC/shifter) on a HPC cluster, as follows:
```
shifter pull jcgoran/coffe:[TAG]
```
Then you can run the code using:
```
shifter run jcgoran/coffe:[TAG] coffe -s [SETTINGSFILE] -n [NUMTHREADS]
```
The above may need some adjustments depending on your HPC environment (using `srun` or similar).

## Testing
To test the output of the code, you can use the command `make check`, which will build the binaries `test_[MODULE]`, where `[MODULE]` can currently be one of `background`, `integrals`, `corrfunc`, `multipoles`, `covariance`, and automatically run them.
Alternatively, you can build them one by one using `make test_[MODULE]`, and run them manually via `./test_[MODULE]`.
This is primarily useful when modifying the code itself, to make sure the old results of the code weren't broken by some new change (feature, bugfix, etc.).

## Bug reports and feature requests
Please use the [issue tracker](https://github.com/JCGoran/coffe/issues) to submit any bug reports and feature requests.
For bug reports, if you are running something other than the Docker version, please specify your platform as well as library versions.

## License and usage conditions
COFFE is licensed under the GNU GPL 3.0. See the `LICENSE` file for more information.
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
