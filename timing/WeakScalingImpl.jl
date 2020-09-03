# Script for processing BATSRUS implicit shocktube test weak scaling runs.
#
# The runlogs are stored in Google Drive, ScalingTests.
#
# Hongyang Zhou, hyzhou@umich.edu 05/09/2019, last modified 09/03/2020

include("Utility.jl")

using PyPlot

Machine = "Frontera"

if Machine == "Bluewaters"
   # New results after modifying the message passing
   result = postprocess("results/gnu_Impl")

   result_0  = sort(result[result.Thread.==0, :],(:Proc))
   result_1  = sort(result[result.Thread.==1, :],(:Proc))
   result_2  = sort(result[result.Thread.==2, :],(:Proc))
   result_4  = sort(result[result.Thread.==4, :],(:Proc))
   result_8  = sort(result[result.Thread.==8, :],(:Proc))
   result_16 = sort(result[result.Thread.==16,:],(:Proc))

   time_0 = copy(result_0.Time)
   time_1 = copy(result_1.Time)
   time_2 = copy(result_2.Time)
   time_4 = copy(result_4.Time)
   time_8 = copy(result_8.Time)
   time_16 = copy(result_16.Time)

   append!(time_0,fill(NaN,size(time_16)[1]-size(time_0)[1]))
   append!(time_1,fill(NaN,size(time_16)[1]-size(time_1)[1]))
   append!(time_2,fill(NaN,size(time_16)[1]-size(time_2)[1]))
   append!(time_4,fill(NaN,size(time_16)[1]-size(time_4)[1]))
   append!(time_8,fill(NaN,size(time_16)[1]-size(time_8)[1]))

   # Calculate the number of cell updates per second for 2nd order scheme
   cell_size = sort(8^3 .* 256 .* result_16.Thread .* result_16.Proc .* 2)

   speed_0 = cell_size ./ time_0
   speed_1 = cell_size ./ time_1
   speed_2 = cell_size ./ time_2
   speed_4 = cell_size ./ time_4
   speed_8 = cell_size ./ time_8
   speed_16 = cell_size ./ time_16

   proc = result_16.Proc .* result_16.Thread
   speed_ref = speed_0[1] ./ proc[1] .* proc

   plt.figure(figsize=(18,5))

   speed = hcat(speed_0, speed_1, speed_2, speed_4, speed_8, speed_16)
   ax1 = plt.subplot(121)
   plt.plot(proc[1:11], speed[1:11,:]/1e9, "-o", markersize=2, alpha=0.8)
   plt.plot(proc[1:11], speed_ref[1:11]/1e9, ":")
   plt.legend(loc="upper left",
         ("pure MPI", "1 thread", "2 threads", "4 threads", "8 threads",
          "16 threads", "linear ref."))
   plt.grid(true)
   ax1.xaxis.set_major_locator(plt.MaxNLocator(5))
   plt.title("BATSRUS Weak Scaling, Implicit Scheme", y=1.03, fontsize=16)
   plt.xlabel("number of cores", fontsize=16)
   plt.ylabel(L"giga cell ($10^9$) updates per second", fontsize=16)

   speed = hcat(speed_0, speed_1, speed_2, speed_4, speed_8, speed_16)
   ax2 = plt.subplot(122)
   plt.plot(proc, speed/1e9, "-o", markersize=3, alpha=0.6)
   plt.plot(proc, speed_ref/1e9, ":")
   plt.legend(loc="upper left",
         ("pure MPI", "1 thread", "2 threads", "4 threads", "8 threads",
          "16 threads", "linear ref."))
   plt.grid(true)
   ax2.xaxis.set_major_locator(plt.MaxNLocator(5))
   plt.title("BATSRUS Weak Scaling, Implicit Scheme", y=1.03, fontsize=16)
   plt.xlabel("number of cores", fontsize=16)
   plt.ylabel(L"giga cell ($10^9$) updates per second", fontsize=16)
   #plt.text(0.95,-0.08, "number of cores", transform=ax2.transAxes)

   #=
   mylabel = ["1 thread", "2 thread", "4 thread", "8 thread", "16 threads",
      "", "", "", "", ""]
   mylabel = reshape(mylabel,1,:)

   p1 = plot(proc,speednew, seriestype=[:line, :scatter],
        markersize=3, lw=1.5, markercolor=[1 2 3 4 5],
        markerstrokewidth = 0.1,
        linecolor=[1 2 3 4 5],
        xscale=:identity, yscale=:identity, legend=:best,
        label=mylabel,
        xlabel="number of cores",
        ylabel="cell updates per second",
        title="BATSRUS Hybrid Weak Scaling for Implicit Scheme",
        xtickfontsize=16, ytickfontsize=16,
        legendfontsize=16, axis=:equal)

   speed_ref = speednew[1] ./ proc[1] .* proc
   plot!(proc,speed_ref, linestyle=:dashdot, linecolor=:black,
         label="linear ref.")
   =#

elseif Machine == "Frontera"
   res = postprocess("RESULT_gnu_Impl")

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

   plt.figure(figsize=(6,5))

   ax1 = plt.subplot(111)
   plt.plot(proc, speed/1e9, "-o", markersize=3, alpha=0.6)
   plt.plot(proc, speed_ref/1e9, ":")
   plt.legend(loc="upper left",
         ("1 thread", "2 threads", "4 threads", "7 threads",
          "14 threads", "28 threads", "56 threads", "linear ref."))
   plt.grid(true)
   ax1.xaxis.set_major_locator(plt.MaxNLocator(5))
   plt.title("BATSRUS Weak Scaling, Implicit Scheme", y=1.03, fontsize=16)
   plt.xlabel("number of cores", fontsize=16)
   plt.ylabel(L"giga cell ($10^9$) updates per second", fontsize=16)
end


