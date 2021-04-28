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

## Guides

This package provides the following functionalities for data analysis in space physics:
  * minimum variance analysis (MVA)
  * spectral analysis
  * moving box average for filtering magnetometer data
  * coordinate transformations


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
