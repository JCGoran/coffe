r"""
## COFFE - COrrelation Function Full-sky Estimator code

This is the Cython (Python) version of the COFFE code, version 3.

### Installation

In brief, if you are on Linux, you can install COFFE using:

```sh
pip install coffe
```

Refer to the installation section of the
[README](https://github.com/JCGoran/coffe/) for installation instructions on
various other platforms.

### Running

#### From a Python environment

For various examples of using COFFE in a Python environment, have a look at the
`coffe.coffe` submodule.

#### From the command line

COFFE also has a command line interface (CLI), which can be invoked after
installation as:

```sh
coffe-cli -s [SETTINGS_FILE] -o [OUTPUT_FILE] -k [TYPE_OF_OUTPUT]
```

where `[SETTINGS_FILE]` is the path to the settings file (for a list of
available settings, as well as defaults, have a look at the [default settings
file](https://github.com/JCGoran/coffe/blob/master/DEFAULT_SETTINGS.cfg)),
`[OUTPUT_FILE]` is the file where the output will be placed, and
`[TYPE_OF_OUTPUT]` is one of `background`, `corrfunc`, `multipoles`, and
`covariance`, for outputs of the background quantities, the 2-point correlation
function, its multipoles, and the covariance of its multipoles, respectively.
"""

from coffe.coffe import Coffe
from coffe.representation import Corrfunc, Covariance, Multipoles
