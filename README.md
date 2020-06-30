# VisAna
[![](https://travis-ci.com/henry2004y/VisAnaJulia.svg?branch=master)][travis-url]
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![](https://img.shields.io/badge/docs-latest-blue.svg)][VisAna-doc]
[![][codecov-img]][codecov-url]

[SWMF](http://csem.engin.umich.edu/tools/swmf/) data reading, converting, visualizing and analyzing using Julia.

For more details, please check the [document][VisAna-doc].

## Prerequisites

Julia 1.4+

## Installation
```
using Pkg
Pkg.add(PackageSpec(url="https://github.com/henry2004y/SWMF", rev="master"))
Pkg.add(PackageSpec(url="https://github.com/henry2004y/VisAnaJulia", rev="master"))
```
Currently Julia does not have a clear way of letting one unregistered package depend on another unregistered package without first add the other package. See the [issue](https://github.com/JuliaLang/Pkg.jl/issues/492) for details. This can be fixed once I register the [SWMF.jl](https://github.com/henry2004y/SWMF.jl) package.

## Usage
```
#using Pkg; Pkg.activate(".") # for dev only
using VisAna
```

See the [examples](docs/src/man/examples.md).

## Guides

This package provides the following functionalities:
  * simulation data reader
  * simulation data visualization
  * data format conversion
  * programming language interoperability
  * data analysis in space physics
  * test particle tracking

The basic functionalities are splitted into a [standalone package](https://github.com/henry2004y/SWMF).
The data analysis part includes spectral analysis, minimum variance analysis and
many functions for aiding the interpretation of data.

See [here](docs/src/man/guide.md) for some development thoughts.
In the future, each part will become a standalone package, and VisAna will only be a container.

## Author

* **Hongyang Zhou** - *Initial work* - [henry2004y](https://github.com/henry2004y)

## Acknowledgments

* All the nice guys who share their codes!

[travis-url]: https://travis-ci.com/henry2004y/VisAnaJulia/builds/
[codecov-img]: https://codecov.io/gh/henry2004y/VisAnaJulia/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/henry2004y/VisAnaJulia
[VisAna-doc]: https://henry2004y.github.io/VisAnaJulia/dev
