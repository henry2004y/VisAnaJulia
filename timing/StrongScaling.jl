# Script for processing the test timing for strong scaling.
#
# Hongyang Zhou, hyzhou@umich.edu 05/09/2019, modified 11/03/2019
include("Utility.jl")

using PyPlot
#using Plots
#pyplot()

# New results after modifying the message passing
res = postprocess("results/StrongScaleRuns/gnu_base32")

res_0 = sort(res[(res.Compiler.=="gnu") .& (res.Thread.==0),:],(:Proc))
res_1 = sort(res[(res.Compiler.=="gnu") .& (res.Thread.==1),:],(:Proc))
res_2 = sort(res[(res.Compiler.=="gnu") .& (res.Thread.==2),:],(:Proc))
res_4 = sort(res[(res.Compiler.=="gnu") .& (res.Thread.==4),:],(:Proc))
res_8 = sort(res[(res.Compiler.=="gnu") .& (res.Thread.==8),:],(:Proc))
res_16 = sort(res[(res.Compiler.=="gnu") .& (res.Thread.==16),:],(:Proc))
res_32 = sort(res[(res.Compiler.=="gnu") .& (res.Thread.==32),:],(:Proc))

# Speed up over 1 node (32 cores) run
speedup_0 = res_0.Time[1] ./ res_0.Time
speedup_1 = res_0.Time[1] ./ res_1.Time
speedup_2 = res_0.Time[1] ./ res_2.Time
speedup_4 = res_0.Time[1] ./ res_4.Time
speedup_8 = res_0.Time[1] ./ res_8.Time
speedup_16 = res_0.Time[1] ./ res_16.Time
speedup_32 = res_0.Time[1] ./ res_32.Time

speedup = hcat(speedup_0, speedup_1, speedup_2, speedup_4, speedup_8,
   speedup_16, speedup_32)

node = res_1.Proc ./ res_1.Proc[1]
proc = res_1.Proc
speed_ref = proc ./ proc[1]


fig = plt.figure(figsize=(5,4))
ax = plt.subplot(1, 1, 1)
plt.plot(proc,speedup,"-o",lw=1,markersize=2,alpha=0.8)
plt.plot(proc,speed_ref,"--")
plt.legend(loc="upper left",
("pure MPI","1 thread", "2 threads", "4 threads", "8 thread", "16 threads",
"32 threads","linear ref."))
plt.grid(true)
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
plt.title("BATSRUS Strong Scaling, Explicit Scheme")
plt.xlabel("number of cores")
plt.ylabel("speedup")
#=
mylabel = ["1 thread", "2 threads", "4 threads", "8 thread", "16 threads",
"32 threads", "Linear ref.", "", "", "", "", "", "", ""]
mylabel = reshape(mylabel, 1, :)

plot(res_1[:Proc]./res_1[:Proc][1],[speedup res_1[:Proc]./res_1[:Proc][1]],
seriestype=[:line, :scatter],
linestyle=[:solid :solid :solid :solid :solid :solid :dash],
markersize=4, lw=2,
legend=:best,
markercolor=[1 2 3 4 5 6 7],
linecolor=[1 2 3 4 5 6 7],
label=mylabel,
xlabel="Number of nodes",
ylabel="Speedup",
title="BATSRUS Hybrid Strong Scaling",
xtickfontsize=16, ytickfontsize=16,
legendfontsize=16,axis=:equal)
#plot!(result[:Proc],result[:Proc])
=#
