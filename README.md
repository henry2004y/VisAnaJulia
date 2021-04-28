# VisAna
[![](https://travis-ci.com/henry2004y/VisAnaJulia.svg?branch=master)][travis-url]
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![](https://img.shields.io/badge/docs-latest-blue.svg)][VisAna-doc]
[![][codecov-img]][codecov-url]

[SWMF](http://csem.engin.umich.edu/tools/swmf/) data reading, converting, visualizing and analyzing using Julia.

For more details, please check the [document][VisAna-doc].

## Prerequisites

Julia 1.6+

## Installation
```
using Pkg
Pkg.add(PackageSpec(url="https://github.com/henry2004y/VisAnaJulia", rev="master"))
```

## Usage

See the [examples](docs/src/man/examples.md).

For local development of the package,
```
using Pkg; Pkg.activate(".")
using VisAna
```

## Guides

This package provides the following functionalities:
  * data analysis in space physics

The basic functionalities are split into a [standalone package](https://github.com/henry2004y/Batsrus.jl).
The data analysis part includes spectral analysis, minimum variance analysis and many functions for aiding the interpretation of data.

See [here](docs/src/man/guide.md) for some development thoughts.
In the future, each part will become a standalone package, and VisAna will only be a container.

## Known Issues

* Currently Julia does not have a clear way of letting one unregistered package depend on another unregistered package without first adding the other package manually. See the [issue](https://github.com/JuliaLang/Pkg.jl/issues/492) for details.

## Author

* **Hongyang Zhou** - *Initial work* - [henry2004y](https://github.com/henry2004y)

## Acknowledgments

* All the nice guys who share their codes!

[travis-url]: https://travis-ci.com/henry2004y/VisAnaJulia/builds/
[codecov-img]: https://codecov.io/gh/henry2004y/VisAnaJulia/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/henry2004y/VisAnaJulia
[VisAna-doc]: https://henry2004y.github.io/VisAnaJulia/dev
