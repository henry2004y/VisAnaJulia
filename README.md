# SpaceAnalysis
[![](https://img.shields.io/github/workflow/status/henry2004y/VisAnaJulia/CI)](https://github.com/henry2004y/VisAnaJulia/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![](https://img.shields.io/badge/docs-latest-blue.svg)][SpaceAnalysis-doc]
[![][codecov-img]][codecov-url]

[SWMF](http://csem.engin.umich.edu/tools/swmf/) data reading, converting, visualizing and analyzing using Julia.

This package provides the following functionalities for data analysis in space physics:
  * minimum variance analysis (MVA)
  * spectral analysis
  * moving box average for filtering magnetometer data
  * coordinate transformations

For more details, please check the [document][SpaceAnalysis-doc].

## Prerequisites

Julia 1.6+

## Installation
```
using Pkg
Pkg.add(PackageSpec(url="https://github.com/henry2004y/VisAnaJulia", rev="master"))
```

## Usage

See the [examples](docs/src/man/examples.md).

## Author

* **Hongyang Zhou** - *Initial work* - [henry2004y](https://github.com/henry2004y)

## Acknowledgments

* All the nice guys who share their codes!

[codecov-img]: https://codecov.io/gh/henry2004y/VisAnaJulia/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/henry2004y/VisAnaJulia
[SpaceAnalysis-doc]: https://henry2004y.github.io/VisAnaJulia/dev
