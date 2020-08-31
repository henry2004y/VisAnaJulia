# Plotting functionalities through Makie.

using Batsrus
using Makie
using AbstractPlotting.MakieLayout

include("utility.jl")

import AbstractPlotting: lines!, contour

function lines!(data::Data, var::AbstractString)
   x, w = data.x, data.w
   VarIndex_ = findindex(data, var)

   # Makie required the same dimension for x and subset of w
   line1 = lines!(ax1, x[:], w[:,VarIndex_], color = :red)
end

function contourf(data::Data, var::AbstractString, levels=0;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1)

   xi, yi, wi = getdata(data, var, plotrange, plotinterval)

   c = contour(xi, yi, wi, fillrange = true)
end

filename = "1d__raw_2_t25.60000_n00000258.out"
data = readdata(filename, dir="../test/data")

#filename = "z=0_raw_1_t25.60000_n00000258.out"
#data = readdata(filename, dir="../test/data")

outer_padding = 30
scene, layout = layoutscene(outer_padding, resolution = (1200, 700),
   backgroundcolor = RGBf0(0.98, 0.98, 0.98))


ax1 = layout[1, 1] = LAxis(scene, title = "p")

lines!(data, "p")


#contourf(data,"p")