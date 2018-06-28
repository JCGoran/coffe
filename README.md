# COFFE
This is the public repository for the code COFFE (COrrelation Function Full-sky Estimator).

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
Extract the contents of the tarball somewhere. Make sure you have the following libraries installed:
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

## Usage
To run the code
```
coffe -s [SETTINGSFILE] -n [NUMTHREADS]
```
where `[SETTINGSFILE]` is the name of the settings file, and `[NUMTHREADS]` is the number of processes you wish to spawn (loosely speaking, the number of cores to use for the computation).

The `settings.cfg` file contains explanations about the possible input and output. For more details, please consult the manual located in the `manual` subdirectory.

## Bug reports and feature requests
Please use the [issue tracker](https://github.com/JCGoran/coffe/issues) to submit any bug reports and feature requests. For bug reports, if you are running something other than the Docker version, please specify your platform and library versions.

## License and usage conditions
COFFE is licensed under the GNU GPL 3.0 license. See the `LICENSE` file for more information.
