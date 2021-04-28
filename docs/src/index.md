# SpaceAnalysis.jl

## Overview

!!! note
    This package is still under development, so be careful for any future breaking changes!

This package provides the following functionalities for data analysis in space physics:
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

## Developers

VisAna is developed by [Hongyang Zhou](https://github.com/henry2004y).
