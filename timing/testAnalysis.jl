# Test script for scalings.
include("Utility.jl")
using Plots

results = postprocess("results/cray_BATL_update/")

compiler= (results[:Compiler].=="cray")

results_1 = sort(results[compiler .& (results[:Thread].==1),:],(:Proc))
results_2 = sort(results[compiler .& (results[:Thread].==2),:],(:Proc))
results_4 = sort(results[compiler .& (results[:Thread].==4),:],(:Proc))
results_8 = sort(results[compiler .& (results[:Thread].==8),:],(:Proc))
results_16 = sort(results[compiler .& (results[:Thread].==16),:],(:Proc))
results_32 = sort(results[compiler .& (results[:Thread].==32),:],(:Proc))

runtime = fill(NaN,13,6)
runtime[1:size(results_1)[1],1] = results_1[:Time]
runtime[1:size(results_2)[1],2] = results_2[:Time]
runtime[1:size(results_4)[1],3] = results_4[:Time]
runtime[1:size(results_8)[1],4] = results_8[:Time]
runtime[1:size(results_16)[1],5] = results_16[:Time]
runtime[1:size(results_32)[1],6] = results_32[:Time]

# Calculate the number of cell updates per second for 2nd order scheme
cell_updates = sort(8^3 .* 256 .* results_32[:Thread] .* results_32[:Proc] .* 2)

speed = hcat([cell_updates ./ runtime[:,i] for i in 1:6]...)

proc = results_32[:Proc] .* results_32[:Thread]


pyplot()
mylabel = reshape(["1 thread", "2 thread", "4 thread", "8 thread", "16 threads",
   "32 threads", "", "", "", "", "", ""],1,12)

#=
plt1 = plot(proc,speed, seriestype=[:line, :scatter], markersize=1, lw=1,
     xscale=:log10, yscale=:log10, legend=:best,
     xlim=[10^1, 10^5.5], ylim=[10^5.5, 10^9.5],
     alpha=0.5,
     markercolor=reshape([(1:13)...],1,13),
     linecolor=reshape([(1:13)...],1,13),
     label=mylabel,
     xlabel="number of cores",
     ylabel="cell updates per second",
     title="BATSRUS Hybrid Weak Scaling",
     xtickfontsize=16, ytickfontsize=16,
     legendfontsize=16)
=#

plt1 = plot(proc,speed, seriestype=[:line, :scatter], markersize=1, lw=1,
          legend=:best,
          markercolor=reshape([(1:13)...],1,13),
          linecolor=reshape([(1:13)...],1,13),
          label=mylabel,
          xlabel="number of cores",
          ylabel="cell updates per second",
          title="BATSRUS Hybrid Weak Scaling",
          xtickfontsize=16, ytickfontsize=16,
          legendfontsize=16)

speed_ref = speed[1] ./ proc[1] .* proc
plot!(proc,speed_ref, linestyle=:dashdot, linecolor=:black,
     label="linear ref.")



results[:core] = results[:Thread] .* results[:Proc]

using StatsPlots
mylabel = reshape(["2^5 cores","2^6 cores","2^7 cores","2^8 cores",
   "2^9 cores","2^10 cores","2^11 cores","2^12 cores","2^13 cores",
   "2^14 cores","2^15 cores","2^16 cores","2^17 cores",
   "","","","","","","","","","","","",""],1,26)

plt2 = plot(reuse = false)
@df results plot(:Thread, :Time, xticks=[2^i for i in 0:5],
   group = (:core),
   label=mylabel,
   seriestype=[:line, :scatter],
   markersize=0.6 .* log2.(:Proc),
   markercolor=reshape([(1:13)...],1,13),
   linecolor=reshape([(1:13)...],1,13),
   xlabel="thread",
   ylabel="execution time [s]",
   title="Timing",
   xtickfontsize=16, ytickfontsize=16,
   legendfontsize=12)
