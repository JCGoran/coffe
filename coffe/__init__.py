r"""
## COFFE - COrrelation Function Full-sky Estimator code

This is the Cython (Python) version of the COFFE code.

### Installation

In brief, if you are on Linux, you can install COFFE using:

```sh
pip install coffe
```

Refer to the [installation
section](https://github.com/JCGoran/coffe/#installation) of the README for
installation instructions on various other platforms.

### Running

For various examples of the capabilities of COFFE, have a look at the
`coffe.coffe` submodule.
"""

from coffe.coffe import Coffe
from coffe.representation import Corrfunc, Covariance, Multipoles
