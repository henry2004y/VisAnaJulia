# VisAna.jl Documentation

## Overview

!!! note
    This package is still under development, so be careful for any future breaking changes!

SWMF data reading, converting, visualizing and analyzing using Julia.

This package is inherited from its predecessor in IDL (developed by G.Tóth) and Matlab (developed by H.Zhou). Currently instead of a real "package", this is more like a collection of scripts. The data loader is split into an stand-alone package [SWMF](https://github.com/henry2004y/SWMF).

This package provides the following functionalities:
  * simulation data reader
  * simulation data visualization
  * data format conversion
  * programming language interoperability
  * data analysis in space physics
  * test particle tracing (WIP)

The data analysis part includes spectral analysis, minimum variance analysis and many functions for aiding the interpretation of data.

The ultimate goal is to build a convenient tool of reading and analyzing simulation outputs which is easy to install, easy to use, and fast.

!!! tip "Ready to use?"
    Feel free to contact the author for any help or collaboration!

## Installation
Install VisAna from the `julia REPL` prompt with
```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/henry2004y/VisAnaJulia", rev="master"))
```

```@contents
Pages = [
    "man/guide.md",
    "man/examples.md",
    "man/functions.md",
    "man/types.md",
    "man/analysis.md"
    "man/testparticle.md"
]
Depth = 1
```

## Benchmark

VisAna has by far the fastest data loading speed among IDL, MATLAB, and Julia.
It has almost the same performance as [spacepy](https://github.com/spacepy/spacepy).
Check the [table](https://github.com/henry2004y/SWMF/blob/master/README.md#Benchmark) for details.

## Developers

VisAna is developed by [Hongyang Zhou](https://github.com/henry2004y).
