# Script for processing BATSRUS explicit shocktube test weak scaling runs.
#
# The runlogs are stored in Google Drive, ScalingTests.
#
# Hongyang Zhou, hyzhou@umich.edu 05/09/2019, last modified 09/03/2020
include("Utility.jl")

using PyPlot
#using Plots
#pyplot()

Machine = "Frontera"

if Machine == "Bluewaters"

   # New results after modifying the message passing
   res = postprocess("results/gnu_BATL_update")
   #res = postprocess("results/gnu")

   result_0 = sort( res[res.Thread.==0,:],(:Proc))
   result_1 = sort( res[res.Thread.==1,:],(:Proc))
   result_2 = sort( res[res.Thread.==2,:],(:Proc))
   result_4 = sort( res[res.Thread.==4,:],(:Proc))
   result_8 = sort( res[res.Thread.==8,:],(:Proc))
   result_16 = sort(res[res.Thread.==16,:],(:Proc))
   result_32 = sort(res[res.Thread.==32,:],(:Proc))

   time_0 = copy(result_0.Time)
   time_1 = copy(result_1.Time)
   time_2 = copy(result_2.Time)
   time_4 = copy(result_4.Time)
   time_8 = copy(result_8.Time)
   time_16 = copy(result_16.Time)
   time_32 = copy(result_32.Time)

   append!(time_0,fill(NaN,size(time_32)[1]-size(time_0)[1]))
   append!(time_1,fill(NaN,size(time_32)[1]-size(time_1)[1]))
   append!(time_2,fill(NaN,size(time_32)[1]-size(time_2)[1]))
   append!(time_4,fill(NaN,size(time_32)[1]-size(time_4)[1]))
   append!(time_8,fill(NaN,size(time_32)[1]-size(time_8)[1]))
   append!(time_16,fill(NaN,size(time_32)[1]-size(time_16)[1]))

   # Calculate the number of cell updates per second for 2nd order scheme
   cell_size = sort(8^3 .* 256 .* result_32.Thread .* result_32.Proc .* 2)

   speed_0 = cell_size ./ time_0
   speed_1 = cell_size ./ time_1
   speed_2 = cell_size ./ time_2
   speed_4 = cell_size ./ time_4
   speed_8 = cell_size ./ time_8
   speed_16 = cell_size ./ time_16
   speed_32 = cell_size ./ time_32

   proc = result_32.Proc .* result_32.Thread

   speed = hcat(speed_0, speed_1, speed_2, speed_4, speed_8, speed_16, speed_32)
   speed_ref = speed_0[1] ./ proc[1] .* proc

   fig = plt.figure(figsize=(18,5))
   # For 1,2,4 threads
   ax1 = plt.subplot(121)

   plt.plot(proc[1:7], speed[1:7,:]/1e9, "-o", lw=1, markersize=3, alpha=0.8)
   plt.plot(proc[1:7], speed_ref[1:7]/1e9, ":")
   plt.legend(loc="upper left",
         ("pure MPI", "1 thread", "2 threads", "4 threads", "8 threads",
         "16 threads", "32 threads", "linear ref."))
   plt.grid(true)
   ax1.xaxis.set_major_locator(plt.MaxNLocator(5))
   plt.title("BATSRUS Weak Scaling, Explicit Scheme", y=1.03, fontsize=16)
   plt.xlabel("number of cores", fontsize=16)
   plt.ylabel(L"giga cell ($10^9$) updates per second", fontsize=16)

   # For 1,2,4,8,16,32 threads
   ax2 = plt.subplot(122)
   plt.plot(proc, speed/1e9, "-o", markersize=3, alpha=0.8)
   plt.plot(proc, speed_ref/1e9, ":")
   plt.legend(loc="upper left",
         ("pure MPI", "1 thread", "2 threads", "4 threads","8 threads",
         "16 threads", "32 threads", "linear ref."))
   plt.grid(true)
   ax2.xaxis.set_major_locator(plt.MaxNLocator(5))
   plt.title("BATSRUS Weak Scaling, Explicit Scheme", y=1.03, fontsize=16)
   plt.xlabel("number of cores", fontsize=16)
   plt.ylabel(L"giga cell ($10^9$) updates per second", fontsize=16)
   #plt.text(0.9,-0.08, "number of cores", transform=ax2.transAxes)


   #=
   mylabel = ["1 thread" "16 threads" "32 threads" "" "" ""]

   plot(proc[1:end],speed[1:end,:], seriestype=[:line, :scatter],
      markersize=4, lw=2, markercolor=[:red :green :purple],
      markerstrokewidth = 0.1,
      linecolor=[:red :green :purple],
      xscale=:identity, yscale=:identity, legend=:best,
      xlim=[10^2, 1.4e5], ylim=[10^6,1.8e9],
      label=mylabel,
      xlabel="number of cores",
      ylabel="cell updates per second",
      title="BATSRUS Hybrid Weak Scaling",
      xtickfontsize=16, ytickfontsize=16,
      legendfontsize=16, axis=:equal,
      size=(900,600))

   plot!(proc[1:end],speed_ref[1:end], linestyle=:dashdot, linecolor=:black,
         label="linear ref.")
   =#

   #=
   speed = hcat(speed_1,speed_2,speed_4,speed_8)
   mylabel = ["1 thread" "2 threads" "4 threads" "8 threads" "" "" "" ""]


   plot(proc[1:6],speed[1:6,:], seriestype=[:line, :scatter],
      markersize=4, lw=2, markercolor=[1 2 3 4],
      markerstrokewidth = 0.1,
      linecolor=[1 2 3 4],
      xscale=:identity, yscale=:identity, legend=:best,
      label=mylabel,
      xlabel="number of cores",
      ylabel="cell updates per second",
      xtickfontsize=16, ytickfontsize=16,
      legendfontsize=16, axis=:equal,
      reuse=false)
   =#

   #using StatsPlots

   #=
   @df results plot(:Thread, :Time, xticks=[1,2,4,8,16,32],
      seriestype=:scatter,
      title="Timing")

   @df results plot!(:Thread, :Time, xticks=[1,2,4,8,16,32],
      seriestype=:line,
      title="Timing")
   =#

   #plotly()

   #=
   @df results plot(:Thread, :Time, xticks=[1,2,4,8,16,32], group=(:Compiler),
      seriestype=[:line,:scatter],
      markersize=0.5 .* log2.(:Proc),
      xlabel="thread",
      ylabel="execution time [s]",
      title="Timing")
   =#

elseif Machine == "Stampede"

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
elseif Machine == "Frontera"

   res = postprocess("RESULT_gnu_Expl")

   result_1 = sort( res[res.Thread.==1,:],(:Proc))
   result_2 = sort( res[res.Thread.==2,:],(:Proc))
   result_4 = sort( res[res.Thread.==4,:],(:Proc))
   result_7 = sort( res[res.Thread.==7,:],(:Proc))
   result_14 = sort( res[res.Thread.==14,:],(:Proc))
   result_28 = sort( res[res.Thread.==28,:],(:Proc))
   result_56 = sort(res[res.Thread.==56,:],(:Proc))

   time_1 = copy(result_1.Time)
   time_2 = copy(result_2.Time)
   time_4 = copy(result_4.Time)
   time_7 = copy(result_7.Time)
   time_14 = copy(result_14.Time)
   time_28 = copy(result_28.Time)
   time_56 = copy(result_56.Time)

   append!(time_1,fill(NaN,size(time_56)[1]-size(time_1)[1]))
   append!(time_2,fill(NaN,size(time_56)[1]-size(time_2)[1]))
   append!(time_4,fill(NaN,size(time_56)[1]-size(time_4)[1]))
   append!(time_7,fill(NaN,size(time_56)[1]-size(time_7)[1]))
   append!(time_14,fill(NaN,size(time_56)[1]-size(time_14)[1]))
   append!(time_28,fill(NaN,size(time_56)[1]-size(time_28)[1]))
   append!(time_56,fill(NaN,size(time_56)[1]-size(time_56)[1]))

   # Calculate the number of cell updates per second for 2nd order scheme
   cell_size = sort(8^3 .* 224 .* result_56.Thread .* result_56.Proc .* 2)

   speed_1 = cell_size ./ time_1
   speed_2 = cell_size ./ time_2
   speed_4 = cell_size ./ time_4
   speed_7 = cell_size ./ time_7
   speed_14 = cell_size ./ time_14
   speed_28 = cell_size ./ time_28
   speed_56 = cell_size ./ time_56

   proc = result_56.Proc .* result_56.Thread

   speed = hcat(speed_1, speed_2, speed_4, speed_7, speed_14, speed_28, speed_56)
   speed_ref = speed_1[1] ./ proc[1] .* proc

   fig = plt.figure(figsize=(6,5))

   ax1 = plt.subplot(111)

   plt.plot(proc, speed/1e9, "-o", lw=1, markersize=3, alpha=0.8)
   plt.plot(proc, speed_ref/1e9, ":")
   plt.legend(loc="upper left",
         ("1 thread", "2 thread", "4 thread", "7 threads", "14 threads", 
         "28 threads", "56 threads", "linear ref."))
   plt.grid(true)
   ax1.xaxis.set_major_locator(plt.MaxNLocator(5))
   plt.title("BATSRUS Weak Scaling, Explicit Scheme", y=1.03, fontsize=16)
   plt.xlabel("number of cores", fontsize=16)
   plt.ylabel(L"giga cell ($10^9$) updates per second", fontsize=16)

end