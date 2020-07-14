# COFFE
This is the public repository for the code COFFE (COrrelation Function Full-sky Estimator).
The relevant theoretical papers are [arXiv:1708.00492](https://arxiv.org/abs/1708.00492) and [arXiv:1806.11090](https://arxiv.org/abs/1806.11090).

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
COFFE has a [Docker](https://docs.docker.com/install/) version available.
To pull the Docker image from DockerHub:
```
docker pull jcgoran/coffe:[TAG]
```
and then run a container using
```
docker run -ti jcgoran/coffe:[TAG]
```
where `[TAG]` is the tag of the release (`latest` by default).
The binary can be run from any directory _inside_ the container.
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
To test the output of the code, you can use use the command `make check`.
Alternatively, you can build tests one by one using `make test_[MODULE]`, where `[MODULE]` can currently be one of `background`, `integrals`, `corrfunc`, `multipoles`, `covariance`.
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
