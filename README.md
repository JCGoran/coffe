# COFFE
This is the public repository for the code COFFE (COrrelation Function Full-sky Estimator).
The relevant theoretical papers are [arXiv:1708.00492](https://arxiv.org/abs/1708.00492) and [arXiv:1806.11090](https://arxiv.org/abs/1806.11090).

## Installation
There are a couple of ways to install the code:

### 1. Docker version
Go [here](https://github.com/JCGoran/coffe/releases) and download the latest release that has "Docker" in the title.
Then install the downloaded image using
```
docker load -i [FILENAME]
```
where `[FILENAME]` is the file you just downloaded. You can then run an interactive session by running
```
docker run -ti coffe-docker:[VERSION]
```
where `[VERSION]` indicates the version of the image you have downloaded. COFFE and all the necessary configuration files are located in the `/build/` directory.

### 2. From the tarball
Go [here](https://github.com/JCGoran/coffe/releases) and download the latest release that doesn't have "Docker" in the title.
Extract the contents of the tarball somewhere. Make sure you have a C compiler (compatible with the C99 standard) and the following libraries installed:
* [FFTW](http://www.fftw.org/download.html)
* [libconfig](https://hyperrealm.github.io/libconfig/)
* GSL
* [CUBA](http://www.feynarts.de/cuba/) (optional)

Once you have installed the dependencies, you can run 
```
./configure && make
```
from the directory in which COFFE was extracted. If you wish to enable the CUBA library, you can specify the flag `--enable-cuba` for the `configure` script.

### 3. From source
After making sure you have the adequate versions of `autotools` installed on your machine, clone this repository and run
```
autoreconf -i
```
From there, you may proceed as under point 2.

### 4. From DockerHub
Run
```
docker pull jcgoran/coffe
```
and then run a container using
```
docker run -ti jcgoran/coffe:latest
```
To generate the binary, run `make`. Copying the input/output can be done using `docker cp` as described in the manual.

## Usage
To run the code
```
coffe -s [SETTINGSFILE] -n [NUMTHREADS]
```
where `[SETTINGSFILE]` is the name of the settings file, and `[NUMTHREADS]` is the number of processes you wish to spawn (loosely speaking, the number of cores to use for the computation).

The `settings.cfg` file contains explanations about the possible input and output. For more details, please consult the manual located in the `manual` subdirectory.

## Bug reports and feature requests
Please use the [issue tracker](https://github.com/JCGoran/coffe/issues) to submit any bug reports and feature requests. For bug reports, if you are running something other than the Docker version, please specify your platform as well as library versions.

## License and usage conditions
COFFE is licensed under the GNU GPL 3.0. See the `LICENSE` file for more information. If you use COFFE in a publication, we kindly ask that you cite the original paper describing the code, located at [arXiv:1806.11090](https://arxiv.org/abs/1806.11090). A `bibTeX` entry is provided below for convenience.
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
