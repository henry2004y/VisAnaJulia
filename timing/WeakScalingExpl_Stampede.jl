# Script for processing BATSRUS explicit shocktube test weak scaling runs.
#
# Hongyang Zhou, hyzhou@umich.edu 05/28/2019
include("Utility.jl")

#using Plots
#pyplot()

using PyPlot

# New results after modifying the message passing
res = postprocess("results/intel_stampede")

result_1 = sort( res[res[:Thread].==1,:],(:Proc))
result_2 = sort( res[res[:Thread].==2,:],(:Proc))
result_4 = sort( res[res[:Thread].==4,:],(:Proc))
result_8 = sort( res[res[:Thread].==8,:],(:Proc))
result_12 = sort( res[res[:Thread].==12,:],(:Proc))
result_16 = sort(res[res[:Thread].==16,:],(:Proc))

time_1 = result_1[:Time]
time_2 = result_2[:Time]
time_4 = result_4[:Time]
time_8 = result_8[:Time]
time_12 = result_12[:Time]
time_16 = result_16[:Time]

append!(time_1,fill(NaN,size(time_16)[1]-size(time_1)[1]))
append!(time_2,fill(NaN,size(time_16)[1]-size(time_2)[1]))
append!(time_4,fill(NaN,size(time_16)[1]-size(time_4)[1]))
append!(time_8,fill(NaN,size(time_16)[1]-size(time_8)[1]))
append!(time_12,fill(NaN,size(time_16)[1]-size(time_12)[1]))
append!(time_16,fill(NaN,size(time_16)[1]-size(time_16)[1]))

# Calculate the number of cell updates per second for 2nd order scheme
cell_size = sort(8^3 .* 256 .* result_16[:Thread] .* result_16[:Proc] .* 2)

speed_1 = cell_size ./ time_1
speed_2 = cell_size ./ time_2
speed_4 = cell_size ./ time_4
speed_8 = cell_size ./ time_8
speed_12 = cell_size ./ time_12
speed_16 = cell_size ./ time_16

proc = result_16[:Proc] .* result_16[:Thread]
#=
speed = hcat(speed_1,speed_16,speed_32)

mylabel = ["1 thread" "16 threads" "32 threads" "" "" ""]

plot(proc,speed, seriestype=[:line, :scatter],
     markersize=4, lw=2, markercolor=[:red :green :purple],
     markerstrokewidth = 0.1,
     linecolor=[:red :green :purple],
     xscale=:identity, yscale=:identity, legend=:best,
     xlim=[10^2, 5.5e5], ylim=[10^6,4.5e9],
     label=mylabel,
     xlabel="number of cores",
     ylabel="cell updates per second",
     title="BATSRUS Hybrid Weak Scaling",
     xtickfontsize=16, ytickfontsize=16,
     legendfontsize=16, axis=:equal)

speed_ref = speed_1[1] ./ proc[1] .* proc
plot!(proc,speed_ref, linestyle=:dashdot, linecolor=:black,
      label="linear ref.")
=#
speed = hcat(speed_1,speed_2,speed_4,speed_8,speed_12,speed_16)
speed_ref = speed_1[1] ./ proc[1] .* proc

plt.figure(figsize=(5,5))
plt.plot(proc,speed/1e9,"-o",markersize=3)
plt.plot(proc,speed_ref/1e9,"--")
plt.legend(loc="upper left",
      ("1 thread", "2 threads", "4 threads", "8 threads", "12 threads",
       "16 threads","linear ref."))
plt.grid(true)
plt.title("BATSRUS Weak Scaling, Explicit Scheme",y=1.03)
plt.xlabel("number of cores")
plt.ylabel(L"giga cell ($10^9$) updates per second")

#=
mylabel = ["1 thread", "2 threads", "4 threads", "8 threads", "16 threads",
      "", "", "", "", ""]
mylabel = reshape(mylabel,1,:)

plot(proc,speed[:,:], seriestype=[:line, :scatter],
     markersize=4, lw=2, markercolor=[1 2 3 4 5],
     markerstrokewidth = 0.1,
     linecolor=[1 2 3 4 5],
     xscale=:identity, yscale=:identity, legend=:best,
     label=mylabel,
     xlabel="number of cores",
     ylabel="cell updates per second",
     xtickfontsize=16, ytickfontsize=16,
     legendfontsize=16, axis=:equal,
     reuse=false)

plot!(proc,speed_ref, linestyle=:dashdot, linecolor=:black,
      label="linear ref.")
=#
