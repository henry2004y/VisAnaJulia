# Examples

Examples of data loading are in [Batsrus.jl](https://henry2004y.github.io/Batsrus.jl/dev/man/examples/).
You can use all the functions in `Batsrus.jl` by, e.g., `VisAna.readdata`, or you can just import the packge by `using Batsrus`.

## Space data analysis

In the [space](../../../space) folder, you can find scripts for comparing magnetic field with observations, cross polar cap potential analysis, diamagnetic current calculation, 1D data frequency analysis, minimum variance analysis, particle phase space distribution plots, cut plots near the X-line reconnection site, and static satellite analysis.

## Plotting with Plots.jl (experimental)

An experimental feature is implemented in [visual_plot.jl](../../../src/visual.jl) for using user recipes in the Julia official plotting package.
It is currently commented out in [VisAna.jl](../../../src/VisAna.jl).
The user recipe allows the plotting functions working on a custom type.
I do not use `Plots.jl` simply because it's too slow and lacks many detailed controls.

## Plotting with Makie (experimental)

Another experimental feature is using Makie for plotting. Makie is known for its GPU support, but the startup time is currently significant, and the default setup is really ugly. However, within years I expect it to become the de facto plotting package. Checkout examples at [visual_Makie.jl](../../../src/visual_Makie.jl).


