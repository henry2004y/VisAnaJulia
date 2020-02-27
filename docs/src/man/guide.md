# Guide

Demos are provided for calling Matlab/Python directly from Julia for debugging and testing. This part will later be separated out for potential Python and Matlab users. Currently the plotting and interpolation needed during plotting are done in Python. For instance, the 3D scatterred interpolation is done via `Interpolate` in Scipy. Hopefully these additional dependencies will be cut down.

The VTK files does not have timestep information. To allow for further time series processing in Paraview, a script `create_pvd.jl` is provided for generating the pvd container.

In principle, I could also try some multi-block (VTM) type for data conversion.

## Tricks

- This is the first time I use Julia for reading general ascii/binary files. It was a pain at first due to the lack of examples and documents using any basic function like read/read!, but fortunately I figured them out myself. One trick in reading binary array data is the usage of view, or subarrays, in Julia. In order to achieve that, I have to implement my own `read!` function in addition to the base ones.
- Tecplot and VTK unstructured data formats have the same connectivity ordering for hexahedron, but different ordering for voxel (in VTK). A function `swaprows` is implemented to switch the orderings.
- Because of the embarrassing parallelism nature of postprocessing, it is quite easy to take advantage of parallel approaches to process the data.
- The built-in streamline function of Matplotlib is not proper for scientifically visualizing field information. The solution is to trace field lines with ODEs and plot the line series, similar to what has been done by [Spacepy](https://github.com/spacepy/spacepy/blob/master/spacepy/pybats/trace2d.py).

## Streamline Tracing

2D streamline tracing is almost finished. 3D streamline tracing is on the way. I need an adaptive step control integration scheme like rk45.

## Particle Tracing

I have a plan of incorporating particle tracing into this module. WIP

## Support for more complicated grid structures

For the plotting, streamline tracing and particle tracing, a common problem is the grid and related interpolation process. I am envisioning a more general approach to deal with block-based and unstructured grid to provide fundamental support for all of these.

##

A real open-source project is a collaborated work not only from a bunch of people, but also a group of languages. In Julia, this can be achieved with the help of the [Package manager](https://julialang.github.io/Pkg.jl/dev/).

I want to have some C dependencies in my code instead of rewriting everything in Julia. This would serve as an attempt to quickly make things work.

Right now this seems to be a little bit difficult for me. I need to learn from experts. The tracing scheme in C is rewritten in Julia so I don't need to bother for now.
Checkout [BinaryBuilder](https://juliapackaging.github.io/BinaryBuilder.jl/latest/#Project-flow-1) for more information. A nice example is given in [this C package](https://github.com/jakubwro/SineWaves.jl).

## Interoperability with Python

In the current version of PyCall and PyJulia, there is already direct support for accessing Julia struct objects (noted as `jlwrap`). 

As to avoid the cross-dependency hail on PyPlot, I split the original package into pure IO [SWMF](https://github.com/henry2004y/SWMF) and post-processing and plotting. This is also a nicer way of organizing larger code base.

## Issues

Before v0.5.1, `readdata` function in Matlab for large data is 2 times faster than that in Julia. The reason is simply using `read` or `unsafe_read` in the lower level. The latter one is much faster. After the fix, Julia version performs 5 times faster than the Matlab version in reading binary data.

At first I forgot to export the Data struct, so everytime when I modified the code and rerun plotdata, it will shout error at me, saying no type was found for the input type.

Precise control of colorbar position in Matplotlib is not an easy task. `axis(“scaled”)` or `axis(“equal”)` will cause issue with the present layout, such as overlapping, cutoff, or too much white spaces. Things are improving, but it takes time. See the scripts in the space folder for some examples of controlling the layouts.

The current support of animation in Matplotlib is not good enough, especially for interactive plotting and scanning through multiple snapshots. The color range is constantly giving me headaches.

The current wrapper over Matplotlib makes it difficult to modify the plots afterwards, which especially causes problems when dealing with time series snapshots. The colorbar is so hard to fix. The solution is, instead of using `level`, provide a range of points.

The support for a long string containing several filenames as inputs has been dropped. It should be substituted by an array of strings.

Right now the derived quantity plots are not supported. In order to achieve this, I may need:
- [x] A new function `get_var(data, filehead, string)` returning the derived variable
- [ ] A new plotting function that understands the derived data type

The first one is achieved by a trick I found on discourse, which basically identifies symbols as names to members in a struct.

There is a user recipe in Plots. This is exactly what I am looking for, but more issues are coming up. I have created a new branch for this development.

I want to do scattered interpolation in Julia directly, but I have not found a simple solution to do this.

A direct wrapper over PyPlot function is possible, and would be more suitable for passing arguments. This may be a more plausible way to go than relying on recipes.

When doing processing in batch mode on a cluster, there's usually no need to render the plots on screen. There exists such a backend for this purpose:
```
using PyPlot
PyPlot.matplotlib.use("Agg")
```
However, notice that currently Agg backend does not support draw_artist. For example, you cannot add an anchored text to your figure.

Vector naming is messed up if you are using Tecplot VTK reader. For example, "B [nT]" --> "B [nT]_X", "B [nT]_Y", "B [nT]_Z". Not a big issue, but annoying.

There is a unit package in Julia [unitful](https://github.com/PainterQubits/Unitful.jl) for handling units. Take a look at that one if you really want to solve the unit problems.

I have encountered a very bad problem of corrupting binary *.vtu files. It turned out that the issue is the starting position of data is wrong because of the way I skip the header AUXDATA part. Sometimes the binary numbers may contain newline character that confuses the reader. It is now fixed. Later on the reading of the header part is completely rewritten to provide better support for a variety of Tecplot Ascii headers.

I have already made a lot of mistakes by mixing the row-major and column-major codes. Explicitly list all the parts that require extra care!

As for the GUI development, GTK seems to be an ideal candidate.

- [x] Fixed colorbar control through Matplotlib
- [x] Test suite for checking validity
- [ ] Full coverage of tests
- [x] Cuts from 3D data visualization besides contour
- [ ] Switch to Makie for 3D plotting and animation
- [ ] PyJulia support for manipulating data directly in Python
- [x] Field tracer 2D in Julia
- [x] Derived variable support
- [x] General postprocessing script for concatenating and converting files.
- [x] Direct wrapper over matplotlib functions to get seamless API
- [x] Replace np.meshgrid with list comprehension
- [ ] Find a substitution of triangulation in Julia
- [ ] Allow dot syntax to get dictionary contents (Base.convert?)
- [ ] Binary library support
- [ ] Macros for quickly looking at data (GUI is the ideal solution!)
- [ ] Magnetic field line plots from simulation
- [ ] Particle phase space distribution plots
