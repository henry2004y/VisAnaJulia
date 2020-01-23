# Examples

IDL format output processing:

- Read data
```
filename = "1d_bin.out";
filehead, data, filelist = readdata(filename);
filehead, data, filelist = readdata(filename, verbose=true);
filehead, data, filelist = readdata(filename, npict=1);
filehead, data, filelist = readdata(filename, dir=".");
```

A general `plotdata` function is provided for quick visualizations. In addition to that, some plotting functions can be directly called as shown below, which allows for more control by the user.

- 1D binary
```
plotdata(data[1], filehead[1], "p", plotmode="line")
plotdata(data[1], filehead[1], "p", plotmode="linegrid")
```

- 2D Cartesian (structured)
```
plotdata(data[1], filehead[1], "p bx;by", plotmode="contbar streamover")
plotdata(data[1], filehead[1], "p bx;by", plotmode="contbar quiverover")
plotdata(data[1], filehead[1], "p bx;by", plotmode="contbar streamover", density=2.0)
plotdata(data[1], filehead[1], "p", plotmode="grid")
plotdata(data[1], filehead[1], "p", plotmode="contbar", plotrange=[-50., 50., -1., 1.])
plotdata(data[1], filehead[1], "p", plotmode="contbar")
plotdata(data[1], filehead[1], "p", plotmode="contbarlog")
plotdata(data[1], filehead[1], "p", plotmode="surfbar")
```

- 2D unstructured
```
plotdata(data[1], filehead[1], "rho", plotmode="contbar")
plotdata(data[1], filehead[1], "rho", plotmode="trimesh")
plotdata(data[1], filehead[1], "rho", plotmode="tricont")
```

- 2D structured spherical coordinates
```
plotdata(data[1], filehead[1], "rho", plotmode="contbar")
```

- 3D box
```
plotdata(data[1], filehead[1], "bx", plotmode="contbar", cut="y", cutPlaneIndex=1, level=20)
plotdata(data[1], filehead[1], "bx", plotmode="contbar", cut="y", plotrange=[-1.4,-1.1,0.70,0.78])
using PyPlot
plt.axis("scaled")
```

- 3D structured spherical coordinates
```
filename = "3d_structured.out";
filehead, data, filelist = readdata(filename, verbose=false);
```

- log file
```
logfilename = "shocktube.log";
filehead, data = readlogdata(logfilename)
```

## Multiple dispatch for matplotlib functions
- line plot
```
plot(data[1], filehead[1], "p", linewidth=2, color="green")
c = plot(data[1], filehead[1], "p")
plt.setp(c, linestyle="--", linewidth=2);
```

- scatter plot
```
scatter(data[1], filehead[1], "p")
```

- contour
```
# 2D contour
contour(data[1], filehead[1], "p")
```

- filled contour
```
# 2D contourf
contourf(data[1], filehead[1], "p")
contourf(data[1], filehead[1], "p", levels, plotrange=[-10,10,-Inf,Inf], plotinterval=0.1)
```

- surface plot
```
# surface
plot_surface(data[1], filehead[1], "p")
```

- triangle surface plot
```
plot_trisurf(data[1], filehead[1], "p")
```

- triangle filled contour plot
```
tricontourf(data[1], filehead[1], "p")
```

- streamline
```
streamplot(data[1], filehead[1], "bx;bz")
streamplot(data[1], filehead[1], "bx;bz", density=2.0, color="k", plotinterval=1.0, plotrange=[-10,10,-Inf,Inf])
```

## Derived variables
```
v = get_vars(data[1], filehead[1], ["Bx", "By", "Bz"])
B = @. sqrt(v.Bx^2 + v.By^2 + v.Bz^2)
```

## Output Format Conversion
ASCII tecplot file:
```
filename = "3d_ascii.dat"
head, data, connectivity  = readtecdata(filename, IsBinary=false)
convertVTK(head, data, connectivity, outname)
```

Binary tecplot file (`DOSAVETECBINARY=TRUE`):
```
filename = "3d_bin.dat"
head, data, connectivity  = readtecdata(filename,true)
convertVTK(head, data, connectivity, outname)
```

Multiple files:
```
using Glob
filenamesIn = "3d*.dat"
dir = "."
filenames = Vector{String}(undef,0)
filesfound = glob(filenamesIn, dir)
filenames = vcat(filenames, filesfound)
tec = readtecdata.(filenames, false) # head, data, connectivity
for (i, outname) in enumerate(filenames)
   convertVTK(tec[i][1], tec[i][2], tec[i][3], outname[1:end-4])
end
```

If each individual file size is large, consider using:
```
using Glob
filenamesIn = "3d*.dat"
dir = "."
filenames = Vector{String}(undef,0)
filesfound = glob(filenamesIn, dir)
filenames = vcat(filenames, filesfound)
for (i, outname) in enumerate(filenames)
   head, data, connectivity = readtecdata(filenames, false)
   convertVTK(head, data, connectivity, outname[1:end-4])
end
```

Multiple files in parallel:
```
using Distributed
@everywhere using Pkg
@everywhere Pkg.activate("VisAnaJulia");
@everywhere using VisAna, Glob

filenamesIn = "cut*.dat"
dir = "."
filenames = Vector{String}(undef,0)
filesfound = glob(filenamesIn, dir)
filenames = vcat(filenames, filesfound)

@sync @distributed for outname in filenames
   println("filename=$(outname)")
   head, data, connectivity = readtecdata(outname, false)
   convertVTK(head, data, connectivity, outname[1:end-4])
end
```