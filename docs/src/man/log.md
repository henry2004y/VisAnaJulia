# Development Log

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
Xi = [y for _ in xi, y in yi]
Yi = [x for x in xi, _ in yi]
wi = interpolator(Xi, Yi)
```

There is a package in Julia called [ScatteredInterpolation.jl](https://github.com/eljungsk/ScatteredInterpolation.jl), but unfortunately it does not have simple bilinear interpolation method, and the existing methods in the package costs too much memory.
From the author of the package:
> I would really like to have an implementation of a fast linear method corresponding to the ones in Python or MATLAB. However, these use a Delaunay triangulation of the sampling points, and as far as I know, the only Julia library providing such a triangulation works only in 2D. I started an implementation two years ago by wrapping the Qhull library that both Python and MATLAB use for the Delaunay triangulation, but that was a major pain and I gave up.

Actually there is a [QHull](https://github.com/JuliaPolyhedra/QHull.jl) wrapper in Julia now. But again this is a wrapper over a Python library. In this case I would say: do not reinvent the wheel for no good reasons.

## Macros

Several places where macros can be used:
* Create a subarray with a name symbol
* Reduce wrapper code duplicates