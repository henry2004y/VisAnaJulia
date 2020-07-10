# Development Log

## Animation

This is a big headache for me right now.
The current support of animation in Matplotlib is not good enough, especially for interactive plotting and scanning through multiple snapshots.

The color range is also an issue.

## Dependency and Package Structure

As to avoid the cross-dependency hail on PyPlot, I split the original package into pure IO [SWMF](https://github.com/henry2004y/SWMF) and post-processing and plotting. This is also a nicer way of organizing larger code base.

Currently VisAna is more of a collection of scripts, instead of a true package.
I am planning to build individual packages for each feature, so that others can make more use of what they want specifically.

Demos are provided for calling Matlab/Python directly from Julia for debugging and testing. This part will later be separated out for potential Python and Matlab users. Currently the plotting and interpolation needed during plotting are done in Python. For instance, the 3D scattered interpolation is done via `Interpolate` in Scipy. Hopefully these additional dependencies will be cut down.

At first I forgot to export the Data struct, so everytime when I modified the code and rerun plotdata, it will shout error at me, saying no type was found for the input type.

Precise control of colorbar position in Matplotlib is not an easy task. `axis(“scaled”)` or `axis(“equal”)` will cause issue with the present layout, such as overlapping, cutoff, or too much white spaces. Things are improving, but it takes time. See the scripts in the space folder for some examples of controlling the layouts.

The current wrapper over Matplotlib makes it difficult to modify the plots afterwards, which especially causes problems when dealing with time series snapshots. The colorbar is so hard to fix. The solution is, instead of using `level`, provide a range of points.

## User Recipe in Plots.jl

There is a user recipe in Plots. This is exactly what I am looking for, but more issues are coming up. I have created a new branch for this development.

I want to do scattered interpolation in Julia directly, but I have not found a simple solution to do this.

## Wrapper over Matplotlib

A direct wrapper over PyPlot function is possible, and would be more suitable for passing arguments. This may be a more plausible way to go than relying on recipes.

When doing processing in batch mode on a cluster, there's usually no need to render the plots on screen. There exists such a backend for this purpose:
```
using PyPlot
PyPlot.matplotlib.use("Agg")
```
However, notice that currently Agg backend does not support draw_artist. For example, you cannot add an anchored text to your figure.

## Streamline

The built-in streamline function of Matplotlib/MATLAB is not proper for scientifically visualizing field information. The solution is to trace field lines with ODEs and plot the line series, similar to what has been done by [Spacepy](https://github.com/spacepy/spacepy/blob/master/spacepy/pybats/trace2d.py).

## GUI

As for the GUI development, GTK seems to be an ideal candidate. However, the interface in Julia lacks full support for the toolkit, which makes it a little bit hard to use.

## Todo List

- [x] Fixed colorbar control through Matplotlib
- [x] Test suite for checking validity
- [ ] Full coverage of tests
- [x] Cuts from 3D data visualization besides contour
- [ ] Switch to Makie for 3D plotting and animation
- [x] Field tracer 2D in Julia
- [x] Derived variable support (dropped because of GUI)
- [x] General postprocessing script for concatenating and converting files.
- [x] Direct wrapper over matplotlib functions to get seamless API
- [x] Replace np.meshgrid with list comprehension
- [ ] Find a substitution of triangulation in Julia
- [ ] Allow dot syntax to get dictionary contents (Base.convert?)
- [ ] Macros for quickly looking at data (GUI is the ideal solution!)
- [x] Magnetic field line plots from simulation
- [x] Particle phase space distribution plots
- [ ] Animation
- [ ] Make more separate small packages instead of one giant collection
