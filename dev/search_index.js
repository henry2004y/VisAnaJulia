var documenterSearchIndex = {"docs":
[{"location":"man/analysis/#Data-Analysis-in-Space-Physics-1","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"","category":"section"},{"location":"man/analysis/#Spectral-Analysis-1","page":"Data Analysis in Space Physics","title":"Spectral Analysis","text":"","category":"section"},{"location":"man/analysis/#FFT-1","page":"Data Analysis in Space Physics","title":"FFT","text":"","category":"section"},{"location":"man/analysis/#Periodogram-1","page":"Data Analysis in Space Physics","title":"Periodogram","text":"","category":"section"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"This is a group of techniques to determine the periodicity of data. Julia has implementations in the DSP package. Here we introduce the usage by looking at practical examples.","category":"page"},{"location":"man/analysis/#Spectrogram-1","page":"Data Analysis in Space Physics","title":"Spectrogram","text":"","category":"section"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"Spectrogram is used a lot in wave analysis. For my purpose, I use it as an approach to visualize time dependent simulation data along a continuous line region.","category":"page"},{"location":"man/analysis/#Minimum-Variance-Analysis-1","page":"Data Analysis in Space Physics","title":"Minimum Variance Analysis","text":"","category":"section"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"A nice introduction is given by Bengt U.Ö.Sonnerup and Maureen Scheible.  Here is a brief summary of the idea. The implementation of MVA can be found in  MVA.jl.","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"The main purpose of minimum or maximum variance analysis (MVA) is to find, from single-spacecraft data, an estimator for the direction normal to a one-dimensional or approximately one-dimensional current layer, wave front, or other transition layer in a plasma.","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"For real transition layers observed in space there are usually more or less pronounced deviations from the ideal 1-D model. The layer is likely to have 2-D or 3-D internal structures which evolve in time and to have temporal fluctuations in the orientation of its normal as well.","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"The minimum variance technique is designed to deal with the situation where some or all of the non-ideal effects mentioned above, except a systematic temporal change in the normal direction, widehatn, are present. As the estimate of  widehatn, the method identifies that direction in space along which the field-component set {mathbfB^(m)cdotwidehatn} (m = 1 2 3M)  has minimum variance. In other words, widehatn is determined by  minimisation of","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"sigma^2 = frac1M sum_m=1^M (mathbfB^(m) - mathbfB)cdotwidehatn ^2","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"where the average langlemathbfBrangle is defined by","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"langlemathbfBrangle equiv frac1M sum_m=1^M mathbfB^(m)","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"and where the minimisation is subject to the normalisation constraint widehatn=1. Using a Lagrange multiplier lambda to implement this constraint, one then seeks the solution of the set of three homogeneous linear equations","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"fracpartialpartial n_xBig( sigma^2 - lambda (widehatn^2 - 1) Big) = 0 \nfracpartialpartial n_yBig( sigma^2 - lambda (widehatn^2 - 1) Big) = 0 \nfracpartialpartial n_zBig( sigma^2 - lambda (widehatn^2 - 1) Big) = 0","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"where sigma^2 is given by the equation above and widehatn is represented in terms of its three components (n_x n_y n_z) along the cartesian coordinate system X, Y, Z (e.g., GSE or GSM) in which the field data mathbfB^(m)  are given. When the differentiations in equations above have been performed, the resulting set of three equations can be written in matrix form as","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"sum_nu=1^3 M_munu^B n_nu = lambda n_mu","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"where the subscripts munu = 123 denote cartesian components along the X, Y, Z system and","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"M_munu^B equiv langle B_mu B_nurangle - langle B_muranglelangle B_nurangle","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"is the magnetic variance matrix. It is seen from the equation that the allowed lambda values are the eigenvalues lambda_1lambda_2lambda_3 (given here in order of decreasing magnitude) of M_munu^B. Since M_munu^B is symmetric, the eigenvalues are all real and the corresponding eigenvectors, x_1, x_2, and x_3, are orthogonal. The three eigenvectors represent the directions of maximum, intermediate, and minimum variance of the field component along each vector.","category":"page"},{"location":"man/analysis/#ULF-Wave-Detection-1","page":"Data Analysis in Space Physics","title":"ULF Wave Detection","text":"","category":"section"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"ULF waves are MHD waves: Alfvén wave, fast wave and slow wave. One basic approach to identify waves is to check the correlation of quantity perturbations.","category":"page"},{"location":"man/analysis/#Correlation-Test-Between-Two-Variables-1","page":"Data Analysis in Space Physics","title":"Correlation Test Between Two Variables","text":"","category":"section"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"This part takes the reference from R.","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"Correlation test is used to evaluate the association between two or more variables.","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"info: Info\n","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"If there is no relationship between the two variables, the average of ``x`` should be the same regardless of ``y`` and vice versa.","category":"page"},{"location":"man/analysis/#Methods-for-correlation-analyses-1","page":"Data Analysis in Space Physics","title":"Methods for correlation analyses","text":"","category":"section"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"There are different methods to perform correlation analysis:","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"Pearson correlation (r), which measures a linear dependence between two variables (x and y). It’s also known as a parametric correlation test because it depends to the distribution of the data. It can be used only when x and y are from normal distribution. The plot of y = f(x) is named the linear regression curve.\nKendall tau and Spearman rho, which are rank-based correlation coefficients (non-parametric).","category":"page"},{"location":"man/analysis/#Correlation-formula-1","page":"Data Analysis in Space Physics","title":"Correlation formula","text":"","category":"section"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"In the formula below,","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"x and y are two vectors of length n\nbarx and bary corresponds to the means of x and y, respectively.","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"Pearson correlation formula","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"r = fracsum (x-barx)(y-barysqrtsum(x-barx)^2sum(y-bary)^2","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"The p-value (significance level) of the correlation can be determined :","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"by using the correlation coefficient table for the degrees of freedom : df=n2, where n is the number of observation in x and y variables.","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"2 or by calculating the t value as follows:","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"t = fracrsqrt1-r^2sqrtn-2","category":"page"},{"location":"man/analysis/#","page":"Data Analysis in Space Physics","title":"Data Analysis in Space Physics","text":"where the corresponding p-value is determined using t table distribution for df=n-2. If the p-value is  5, then the correlation between x and y is significant.","category":"page"},{"location":"man/guide/#Guide-1","page":"Guide","title":"Guide","text":"","category":"section"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"Demos are provided for calling Matlab/Python directly from Julia for debugging and testing. This part will later be separated out for potential Python and Matlab users. Currently the plotting and interpolation needed during plotting are done in Python. For instance, the 3D scatterred interpolation is done via Interpolate in Scipy. Hopefully these additional dependencies will be cut down.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"The VTK files does not have timestep information. To allow for further time series processing in Paraview, a script create_pvd.jl is provided for generating the pvd container.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"In principle, I could also try some multi-block (VTM) type for conversion.","category":"page"},{"location":"man/guide/#Tricks-1","page":"Guide","title":"Tricks","text":"","category":"section"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"This is the first time I use Julia for reading general ascii/binary files. It was a pain at first due to the lack of examples and documents using any basic function like read/read!, but fortunately I figured them out myself. One trick in reading binary array data is the usage of view, or subarrays, in Julia. In order to achieve that, I have to implement my own read! function in addition to the base ones.\nTecplot and VTK unstructured data formats have the same connectivity ordering for hexahedron, but different ordering for voxel (in VTK). A function swaprows is implemented to switch the orderings.\nBecause of the embarrassing parallelism nature of postprocessing, it is quite easy to take advantage of parallel approaches to process the data.","category":"page"},{"location":"man/guide/#Issues-1","page":"Guide","title":"Issues","text":"","category":"section"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"At first I forgot to export the Data struct, so everytime when I modified the code and rerun plotdata, it will shout error at me, saying no type was found for the input type.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"The current support of animation in Matplotlib is not good enough, especially for interactive plotting and scanning through multiple snapshots. The color range is constantly giving me headaches.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"The current wrapper over Matplotlib makes it difficult to modify the plots afterwards, which especially causes problems when dealing with time series snapshots. The colorbar is so hard to fix.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"In the roadmap of PyCall 2.0, there will direct support for accessing Julia objects. I hesitate to do it myself, so let's just wait for it to come.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"The support for a long string containing several filenames as inputs has been dropped. It should be substituted by an array of strings.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"Right now the derived quantity plots are not supported. In order to achieve this, I may need:","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"[x] A new function get_var(data, filehead, string) returning the derived variable\n[ ] A new plotting function that understands the derived data type","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"The first one is achieved by a trick I found on discourse, which basically identifies symbols as names to members in a struct.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"There is a user recipe in Plots. This is exactly what I am looking for, but more issues are coming up. I have created a new branch for this development.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"I want to do scattered interpolation in Julia directly, but I have not found a simple solution to do this.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"A direct wrapper over PyPlot function is possible, and would be more suitable for passing arguments. This may be a more plausible way to go than relying on recipes.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"When doing processing in batch mode on a cluster, there's usually no need to render the plots on screen. There exists such a backend for this purpose:","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"using PyPlot\nPyPlot.matplotlib.use(\"Agg\")","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"However, notice that currently Agg backend does not support draw_artist. For example, you cannot add an anchored text to your figure.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"Vector naming is messed up if you are using Tecplot VTK reader. For example, \"B [nT]\" –> \"B [nT]X\", \"B [nT]Y\", \"B [nT]_Z\". Not a big issue, but annoying.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"There is a unit package in Julia unitful for handling units. Take a look at that one if you really want to solve the unit problems.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"I have encountered a very bad problem of corrupting binary *.vtu files. It turned out that the issue is the starting position of data is wrong because of the way I skip the header AUXDATA part. Sometimes the binary numbers may contain newline character that confuses the reader. It is now fixed.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"[x] Fixed colorbar control through Matplotlib\n[x] Test suite for checking validity\n[ ] Full coverage of tests\n[x] Cuts from 3D data visualization besides contour\n[ ] Switch to Makie for 3D plotting and animation\n[ ] PyBase support for manipulating data directly in Python\n[x] Derived variable support\n[ ] General postprocessing script for concatenating and converting files.\n[x] Direct wrapper over matplotlib functions to get seamless API\n[x] Replace np.meshgrid with list comprehension\n[ ] Find a substitution of triangulation in Julia\n[ ] Allow dot syntax to get dictionary contents (Base.convert?)","category":"page"},{"location":"man/examples/#Examples-1","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"IDL format output processing:","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"Read data","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"filename = \"1d_bin.out\";\nfilehead, data, filelist = readdata(filename);\nfilehead, data, filelist = readdata(filename, verbose=true);\nfilehead, data, filelist = readdata(filename, npict=1);\nfilehead, data, filelist = readdata(filename, dir=\".\");","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"A general plotdata function is provided for quick visualizations. In addition to that, some plotting functions can be directly called as shown below, which allows for more control by the user.","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"1D binary","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"plotdata(data[1], filehead[1], \"p\", plotmode=\"line\")\nplotdata(data[1], filehead[1], \"p\", plotmode=\"linegrid\")","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"2D Cartesian (structured)","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"plotdata(data[1], filehead[1], \"p bx;by\", plotmode=\"contbar streamover\")\nplotdata(data[1], filehead[1], \"p bx;by\", plotmode=\"contbar quiverover\")\nplotdata(data[1], filehead[1], \"p bx;by\", plotmode=\"contbar streamover\", density=2.0)\nplotdata(data[1], filehead[1], \"p\", plotmode=\"grid\")\nplotdata(data[1], filehead[1], \"p\", plotmode=\"contbar\", plotrange=[-50., 50., -1., 1.])\nplotdata(data[1], filehead[1], \"p\", plotmode=\"contbar\")\nplotdata(data[1], filehead[1], \"p\", plotmode=\"contbarlog\")\nplotdata(data[1], filehead[1], \"p\", plotmode=\"surfbar\")","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"2D unstructured","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"plotdata(data[1], filehead[1], \"rho\", plotmode=\"contbar\")\nplotdata(data[1], filehead[1], \"rho\", plotmode=\"trimesh\")\nplotdata(data[1], filehead[1], \"rho\", plotmode=\"tricont\")","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"2D structured spherical coordinates","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"plotdata(data[1], filehead[1], \"rho\", plotmode=\"contbar\")","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"3D box","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"plotdata(data[1], filehead[1], \"bx\", plotmode=\"contbar\", cut=\"y\", cutPlaneIndex=1, level=20)\nplotdata(data[1], filehead[1], \"bx\", plotmode=\"contbar\", cut=\"y\", plotrange=[-1.4,-1.1,0.70,0.78])\nusing PyPlot\nplt.axis(\"scaled\")","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"3D structured spherical coordinates","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"filename = \"3d_structured.out\";\nfilehead, data, filelist = readdata(filename, verbose=false);","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"log file","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"logfilename = \"shocktube.log\";\nfilehead, data = readlogdata(logfilename)","category":"page"},{"location":"man/examples/#Multiple-dispatch-for-matplotlib-functions-1","page":"Examples","title":"Multiple dispatch for matplotlib functions","text":"","category":"section"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"line plot","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"plot(data[1], filehead[1], \"p\", linewidth=2, color=\"green\")\nc = plot(data[1], filehead[1], \"p\")\nplt.setp(c, linestyle=\"--\", linewidth=2);","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"scatter plot","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"scatter(data[1], filehead[1], \"p\")","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"contour","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"# 2D contour\ncontour(data[1], filehead[1], \"p\")","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"filled contour","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"# 2D contourf\ncontourf(data[1], filehead[1], \"p\")\ncontourf(data[1], filehead[1], \"p\", levels, plotrange=[-10,10,-Inf,Inf], plotinterval=0.1)","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"surface plot","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"# surface\nplot_surface(data[1], filehead[1], \"p\")","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"triangle surface plot","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"plot_trisurf(data[1], filehead[1], \"p\")","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"triangle filled contour plot","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"tricontourf(data[1], filehead[1], \"p\")","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"streamline","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"streamplot(data[1], filehead[1], \"bx;bz\")\nstreamplot(data[1], filehead[1], \"bx;bz\", density=2.0, color=\"k\", plotinterval=1.0, plotrange=[-10,10,-Inf,Inf])","category":"page"},{"location":"man/examples/#Derived-variables-1","page":"Examples","title":"Derived variables","text":"","category":"section"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"v = get_vars(data[1], filehead[1], [\"Bx\", \"By\", \"Bz\"])\nB = @. sqrt(v.Bx^2 + v.By^2 + v.Bz^2)","category":"page"},{"location":"man/examples/#Output-Format-Conversion-1","page":"Examples","title":"Output Format Conversion","text":"","category":"section"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"ASCII tecplot file:","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"filename = \"3d_ascii.dat\"\nhead, data, connectivity  = readtecdata(filename,false)\nconvertVTK(head, data, connectivity, outname)","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"Binary tecplot file (DOSAVETECBINARY=TRUE):","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"filename = \"3d_bin.dat\"\nhead, data, connectivity  = readtecdata(filename,true)\nconvertVTK(head, data, connectivity, outname)","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"Multiple files:","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"using Glob\nfilenamesIn = \"3d*.dat\"\ndir = \".\"\nfilenames = Vector{String}(undef,0)\nfilesfound = glob(filenamesIn, dir)\nfilenames = vcat(filenames, filesfound)\ntec = readtecdata.(filenames, false) # head, data, connectivity\nfor (i, outname) in enumerate(filenames)\n   convertVTK(tec[i][1], tec[i][2], tec[i][3], outname[1:end-4])\nend","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"If each individual file size is large, consider using:","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"using Glob\nfilenamesIn = \"3d*.dat\"\ndir = \".\"\nfilenames = Vector{String}(undef,0)\nfilesfound = glob(filenamesIn, dir)\nfilenames = vcat(filenames, filesfound)\nfor (i, outname) in enumerate(filenames)\n   head, data, connectivity = readtecdata(filenames, false)\n   convertVTK(head, data, connectivity, outname[1:end-4])\nend","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"Multiple files in parallel:","category":"page"},{"location":"man/examples/#","page":"Examples","title":"Examples","text":"using Distributed\n@everywhere using Pkg\n@everywhere Pkg.activate(\"VisAnaJulia\");\n@everywhere using VisAna, Glob\n\nfilenamesIn = \"cut*.dat\"\ndir = \".\"\nfilenames = Vector{String}(undef,0)\nfilesfound = glob(filenamesIn, dir)\nfilenames = vcat(filenames, filesfound)\n\n@sync @distributed for outname in filenames\n   println(\"filename=$(outname)\")\n   head, data, connectivity = readtecdata(outname, false)\n   convertVTK(head, data, connectivity, outname[1:end-4])\nend","category":"page"},{"location":"#VisAna.jl-Documentation-1","page":"VisAna.jl Documentation","title":"VisAna.jl Documentation","text":"","category":"section"},{"location":"#Overview-1","page":"VisAna.jl Documentation","title":"Overview","text":"","category":"section"},{"location":"#","page":"VisAna.jl Documentation","title":"VisAna.jl Documentation","text":"note: Note\nThis package is still under development, so be careful for any future breaking changes!","category":"page"},{"location":"#","page":"VisAna.jl Documentation","title":"VisAna.jl Documentation","text":"SWMF data reader and visualization using Julia.","category":"page"},{"location":"#","page":"VisAna.jl Documentation","title":"VisAna.jl Documentation","text":"This package is the inherited from its predecessor in IDL (developed by G.Toth) and Matlab (developed by H.Zhou). It can be combined with the VTK format converter writeVTK to generate files for Paraview and Tecplot. By default the file size will be reduced with compression level 6, but the actual compression ratio depends on the original data.","category":"page"},{"location":"#","page":"VisAna.jl Documentation","title":"VisAna.jl Documentation","text":"This package consists of five parts:","category":"page"},{"location":"#","page":"VisAna.jl Documentation","title":"VisAna.jl Documentation","text":"simulation data reader\nsimulation data visualization\ndata format conversion\nprogramming language interoperability\ndata analysis in space physics","category":"page"},{"location":"#","page":"VisAna.jl Documentation","title":"VisAna.jl Documentation","text":"The data analysis part includes spectral analysis, minimum variance analysis and many functions for aiding the interpretation of data.","category":"page"},{"location":"#","page":"VisAna.jl Documentation","title":"VisAna.jl Documentation","text":"The ultimate goal is to build a convenient tool of reading and analyzing simulation outputs which is easy to install and easy to use.","category":"page"},{"location":"#","page":"VisAna.jl Documentation","title":"VisAna.jl Documentation","text":"tip: Ready to use?\nFeel free to contact the author for any help or collaboration!","category":"page"},{"location":"#Installation-1","page":"VisAna.jl Documentation","title":"Installation","text":"","category":"section"},{"location":"#","page":"VisAna.jl Documentation","title":"VisAna.jl Documentation","text":"Install VisAna from the julia REPL prompt with","category":"page"},{"location":"#","page":"VisAna.jl Documentation","title":"VisAna.jl Documentation","text":"using Pkg\nPkg.add(PackageSpec(url=\"https://github.com/henry2004y/VisAnaJulia\", rev=\"master\"))","category":"page"},{"location":"#","page":"VisAna.jl Documentation","title":"VisAna.jl Documentation","text":"Pages = [\n    \"man/guide.md\",\n    \"man/examples.md\",\n    \"man/functions.md\",\n    \"man/types.md\",\n    \"man/analysis.md\"\n]\nDepth = 1","category":"page"},{"location":"#Developers-1","page":"VisAna.jl Documentation","title":"Developers","text":"","category":"section"},{"location":"#","page":"VisAna.jl Documentation","title":"VisAna.jl Documentation","text":"VisAna is developed by Hongyang Zhou.","category":"page"},{"location":"man/types/#Private-types-1","page":"Private types","title":"Private types","text":"","category":"section"},{"location":"man/types/#Private-types-in-module-VisAna:-1","page":"Private types","title":"Private types in module VisAna:","text":"","category":"section"},{"location":"man/types/#","page":"Private types","title":"Private types","text":"Modules = [VisAna]\nPublic = false\nOrder = [:type]","category":"page"},{"location":"man/functions/#Functions-1","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"man/functions/#Functions-exported-from-VisAna:-1","page":"Functions","title":"Functions exported from VisAna:","text":"","category":"section"},{"location":"man/functions/#","page":"Functions","title":"Functions","text":"Modules = [VisAna]\nPrivate = false\nOrder = [:function]","category":"page"},{"location":"man/functions/#PyPlot.contour","page":"Functions","title":"PyPlot.contour","text":"contour(data, filehead, var, levels=0; plotrange, plotinterval, kwargs)\n\nWrapper over the contour function in matplotlib.\n\n\n\n\n\n","category":"function"},{"location":"man/functions/#PyPlot.contourf","page":"Functions","title":"PyPlot.contourf","text":"contourf(data, filehead, var, levels=0; plotrange, plotinterval, kwargs)\n\nWrapper over the contourf function in matplotlib.\n\n\n\n\n\n","category":"function"},{"location":"man/functions/#PyPlot.plot-Tuple{Data,Dict,String}","page":"Functions","title":"PyPlot.plot","text":"plot(data, filehead, var; kwargs)\n\nWrapper over the plot function in matplotlib.\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#PyPlot.plot_surface-Tuple{Data,Dict,String}","page":"Functions","title":"PyPlot.plot_surface","text":"plot_surface(data, filehead, var; plotrange, plotinterval, kwargs)\n\nWrapper over the plot_surface function in matplotlib.\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#PyPlot.plot_trisurf-Tuple{Data,Dict,String}","page":"Functions","title":"PyPlot.plot_trisurf","text":"plot_trisurf(data::Data, filehead::Dict, var::String;\n\tplotrange::Vector{Float64}=[-Inf,Inf,-Inf,Inf], kwargs::Dict=Dict())\n\nWrapper over the plot_trisurf function in matplotlib.\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#PyPlot.scatter-Tuple{Data,Dict,String}","page":"Functions","title":"PyPlot.scatter","text":"scatter(data, filehead, var; kwargs)\n\nWrapper over the scatter function in matplotlib.\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#PyPlot.streamplot-Tuple{Data,Dict,String}","page":"Functions","title":"PyPlot.streamplot","text":"streamplot(data, filehead, var; plotrange, plotinterval)\n\nWrapper over the streamplot function in matplotlib. Streamplot does not have **kwargs in the API.\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#PyPlot.tricontourf-Tuple{Data,Dict,String}","page":"Functions","title":"PyPlot.tricontourf","text":"tricontourf(data, filehead, var; plotrange, plotinterval, kwargs)\n\nWrapper over the tricontourf function in matplotlib.\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#VisAna.animatedata-Tuple{FileList,String}","page":"Functions","title":"VisAna.animatedata","text":"animatedata(filelist, func, (plotmode=\"contbar\",\n  plotrange=[-Inf Inf -Inf Inf],\n  plotinterval=0.1))\n\nGenerate animations from data. This is basically calling plotdata function for multiple snapshots. The main issue here is to determine the colorbar/axis range in advance to avoid any jump in the movie.\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#VisAna.convertVTK","page":"Functions","title":"VisAna.convertVTK","text":"convertVTK(head, data, connectivity, filename)\n\nConvert 3D unstructured Tecplot data to VTK. Note that if using voxel type data in VTK, the connectivity sequence is different from Tecplot.\n\n\n\n\n\n","category":"function"},{"location":"man/functions/#VisAna.plotdata-Tuple{Data,Dict,String}","page":"Functions","title":"VisAna.plotdata","text":"plotdata(data, filehead, func, (...))\n\nPlot the variable from SWMF output.\n\nplotdata(data, filehead, \"p\", plotmode=\"contbar\")\n\nplotdata(data, filehead, \"p\", plotmode=\"grid\")\n\nplotdata(data, filehead, func, plotmode=\"trimesh\",plotrange=plotrange,    plotinterval=0.2)\n\nInput arguments\n\ndata::Data: original variable data.\nfilehead::Dict: header information.\nvars::String: variables for plotting.\nplotmode::String: (optional) type of plotting [\"cont\",\"contbar\"]...\nplotrange::Vector: (optional) range of plotting.\nplotinterval: (optional) interval for interpolation.\nlevel: (optional) level of contour.\ndensity: (optional) density for streamlines.\ncut: (optional) select 2D cut plane from 3D outputs [\"x\",\"y\",\"z\"].\ncutPlaneIndex: (optional)\nmultifigure: (optional) 1 for multifigure display, 0 for subplots.\nverbose: (optional) display additional information.\n\nRight now this can only deal with 2D plots or 3D cuts. Full 3D plots may be supported in the future.\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#VisAna.plotlogdata-Tuple{Data,Dict,String}","page":"Functions","title":"VisAna.plotlogdata","text":"plotlogdata(data, filehead, vars, (plotmode=\"line\", plotrange=[-Inf,Inf]))\n\nPlot information from log file.\n\nInput arguments\n\ndata::Data: original variable data.\nfilehead::Dict: header information.\nvars::String: variables for plotting.\nplotmode::String: (optional) type of plotting [\"line\",\"scatter\"].\nplotrange::Vector: (optional) range of plotting.\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#VisAna.readdata-Tuple{String}","page":"Functions","title":"VisAna.readdata","text":"readdata(filenames,(, dir=\".\", npict=1, verbose=false))\n\nRead data from BATSRUS output files. Stores the npict-th snapshot from an ascii or binary data file into the x [coordinates] and w [data] arrays. Filenames can be provided with wildcards.\n\nExamples\n\nfilenames = \"1d__raw*\"\nfileheads, data, filelist = readdata(filenames)\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#VisAna.readlogdata-Tuple{String}","page":"Functions","title":"VisAna.readlogdata","text":"readlogdata(filename)\n\nRead information from log file.\n\n\n\n\n\n","category":"method"},{"location":"man/functions/#VisAna.readtecdata","page":"Functions","title":"VisAna.readtecdata","text":"readtecdata(filename, IsBinary=false, verbose=false)\n\nReturn header, data and connectivity from BATSRUS Tecplot outputs. Both binary and ascii formats are supported.\n\nExamples\n\nfilename = \"3d_ascii.dat\"\nhead, data, connectivity = readtecdata(filename)\n\n\n\n\n\n","category":"function"}]
}
