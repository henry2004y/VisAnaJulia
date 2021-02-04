# VisAna.jl Documentation

## Overview

!!! note
    This package is still under development, so be careful for any future breaking changes!

[SWMF](http://csem.engin.umich.edu/tools/swmf) data reading, converting, visualizing and analyzing using Julia.

This package is inherited from its predecessor in IDL (developed by G.Tóth) and Matlab (developed by H.Zhou). Currently instead of a real "package", this is more like a collection of scripts. The data loader and converter is split into an stand-alone package [Batsrus.jl](https://github.com/henry2004y/Batsrus.jl).

This package provides the following functionalities:
  * SWMF simulation data visualization
    * phase space distribution
    * 2D slices of 3D data, including all common plots
    * log variable plots
  * data analysis in space physics
    * magnetic field comparison with satellite data
    * spectral analysis
    * minimum variance analysis
    * wave analysis

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
    "man/log.md",
    "man/examples.md",
    "man/functions.md",
    "man/types.md",
    "man/analysis.md"
    "man/trace.md"
    "man/text.md"
]
Depth = 1
```

## Benchmark

VisAna has by far the fastest data loading speed among IDL, MATLAB, and Julia.
It has almost the same performance as [spacepy](https://github.com/spacepy/spacepy).
Check the [table](https://henry2004y.github.io/Batsrus.jl/dev/#Benchmark) for details.

## Developers

VisAna is developed by [Hongyang Zhou](https://github.com/henry2004y).
