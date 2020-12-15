# Development Log

## Animation

This is a big headache for me right now.
The current support of animation in Matplotlib is not good enough, especially for interactive plotting and scanning through multiple snapshots.

The color range is also an issue.

## Dependency and Package Structure

As to avoid the cross-dependency hell on PyPlot, I split the original package into pure IO [Batsrus.jl](https://github.com/henry2004y/Batsrus.jl) and post-processing and plotting. This is also a nicer way of organizing larger code base.
Furthermore, after Julia 1.6, the Python style `import A as B` and `using LinearAlgebra: cholesky as c, lu as l` are supported.
This makes it possible to completely drop the dependency on PyPlot.jl and instead everytime one would need to say `using VisAna, PyPlot` to trigger the Matplotlib functions. A very nice thing is that I can in principle switch between PyPlot.jl and Plots.jl easily, depending on which package to use!

Currently VisAna is more of a collection of scripts, instead of a true package.
I am planning to build individual packages for each feature, so that others can make more use of what they want specifically.

Demos are provided for calling Matlab/Python directly from Julia for debugging and testing. This part will later be separated out for potential Python and Matlab users. Currently the plotting and interpolation needed during plotting are done in Python. For instance, the scattered interpolation is done via `Interpolate` in Scipy. Hopefully these additional dependencies will be cut down.

At first I forgot to export the Data struct, so everytime when I modified the code and rerun plotdata, it will shout error at me, saying no type was found for the input type.

Precise control of colorbar position in Matplotlib is not an easy task. `axis(“scaled”)` or `axis(“equal”)` will cause issue with the present layout, such as overlapping, cutoff, or too much white spaces. Things are improving, but it takes time. See the scripts in the space folder for some examples of controlling the layouts.

The current wrapper over Matplotlib makes it difficult to modify the plots afterwards, which especially causes problems when dealing with time series snapshots. The colorbar is so hard to fix. The solution is, instead of using `level`, provide a range of points.

## User Recipe in Plots.jl

There is a *extremely powerful* user recipe in Plots.

* Repeatly using the same GKSTerm on Mac will display only white backgrounds in the end.
* By default Plots uses `gr()` backend. The GR backend contour plot only accept vector x,y!
* I don't want to have `Plots.jl` as a dependency. With simple plotting features this can work, but we may encounter [issues](https://github.com/JuliaPlots/RecipesBase.jl/issues/72) later. After Julia 1.6 this may be completely solved!
* There is already a [UnitfulRecipes.jl](https://github.com/jw3126/UnitfulRecipes.jl) that provides the capability of auto-displaying units in figure labels, and it works smoothly with my user recipe. Amazing.
* I have already built a customized package [UnitfulBatsrus.jl](https://github.com/henry2004y/UnitfulBatsrus.jl.git) and set it as a dependency for VisAna. Instead of the usual `u"km/s"` notation, we just need to use `bu"amucc"`.


## Scattered interpolation

SWMF outputs may be in generalized coordinates. For the purpose of plotting, we often need to first interpolate onto a uniform mesh.
Even thought this seems to be a very basic job, I have not yet found a simple solution to do this in native Julia.
My current workaround is to call Python:
```
using PyPlot
n = 100
X, Y, W = rand(n), rand(n), rand(n)
interval = 0.02
xi = range(minimum(X), stop=maximum(X), step=interval)
yi = range(minimum(Y), stop=maximum(Y), step=interval)
# Perform linear interpolation of the data (x,y) on grid(xi,yi)
triang = matplotlib.tri.Triangulation(X,Y)
interpolator = matplotlib.tri.LinearTriInterpolator(triang, W)
Xi = [y for x in xi, y in yi]
Yi = [x for x in xi, y in yi]
wi = interpolator(Xi, Yi)
```

There is a package in Julia called [ScatteredInterpolation.jl](https://github.com/eljungsk/ScatteredInterpolation.jl), but unfortunately it does not have simple bilinear interpolation method, and the existing methods in the package costs too much memory.
From the author of the package:
> I would really like to have an implementation of a fast linear method corresponding to the ones in Python or MATLAB. However, these use a Delaunay triangulation of the sampling points, and as far as I know, the only Julia library providing such a triangulation works only in 2D. I started an implementation two years ago by wrapping the Qhull library that both Python and MATLAB use for the Delaunay triangulation, but that was a major pain and I gave up.

Actually there is a [QHull](https://github.com/JuliaPolyhedra/QHull.jl) wrapper in Julia now. But again this is a wrapper over a Python library. In this case I would say: do not reinvent the wheel for no good reasons.

## Wrapper over Matplotlib

A direct wrapper over PyPlot function is possible, and would be more suitable for passing arguments. This may be a more plausible way to go than relying on recipes.

When doing processing in batch mode on a cluster, there's usually no need to render the plots on screen. There exists such a backend for this purpose:
```
using PyPlot
PyPlot.matplotlib.use("Agg")
```
However, notice that currently Agg backend does not support draw_artist. For example, you cannot add an anchored text to your figure.

Unlike the user recipes in `Plots.jl`, using `PyPlot.jl` would require to have it as a dependency. (This might not be necessary for the upcoming Julia 1.6!)

## Makie

[MakieLayout](https://jkrumbiegel.github.io/MakieLayout.jl/dev/) is a nice extension built on top of Makie to create publication quality figures and interactive plots.
It basically includes all the funcationalities I want, so definitely worth a try.

## Macros

Several places where macros can be used:
* Create a subarray with a name symbol
* Reduce wrapper code duplicates

## Streamline

The built-in streamline function of Matplotlib/MATLAB is not proper for scientifically visualizing field information. The solution is to trace field lines with ODEs and plot the line series, similar to what has been done by [Spacepy](https://github.com/spacepy/spacepy/blob/master/spacepy/pybats/trace2d.py).

## GUI

As for the GUI development, GTK seems to be an ideal candidate. However, the [GTK interface in Julia](https://github.com/JuliaGraphics/Gtk.jl) lacks full support for the toolkit, which makes it a little bit hard to use. I have only played with it for half a day. You can design the appearance of your window interactively, and save your in an HTML-like file.

At this point GUI is not necessarily needed, if it does not speed up my own workflow.

## Todo List

- [x] Fixed colorbar control through Matplotlib
- [x] Test suite for checking validity
- [ ] Full coverage of tests
- [x] Cuts from 3D data visualization besides contour
- [ ] Port to Makie
- [x] Field tracer 2D in Julia
- [x] Derived variable support
- [x] General postprocessing script for concatenating and converting files.
- [x] Direct wrapper over Matplotlib functions to get seamless API
- [x] Replace np.meshgrid with list comprehension
- [ ] Find a substitution of triangulation in Julia
- [ ] Allow dot syntax to get dictionary contents (Base.convert?)
- [ ] Macros for quickly looking at data (GUI is the ideal solution!)
- [x] Magnetic field line plots from simulation
- [x] Particle phase space distribution plots
- [ ] Animation
- [ ] Make more separate small packages instead of one giant collection
