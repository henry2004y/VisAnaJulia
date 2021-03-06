# Script for SWMF Galileo flybys simulated magnetic field comparison.
#
# 1. Steady state
# 2. Time accurate
#
# Hongyang Zhou, hyzhou@umich.edu 02/03/2020

using Batsrus, CSV, Dates, PyPlot, Interpolations

struct BField{T<:AbstractFloat}
   t::Vector{DateTime}
   Bx::AbstractVector{T}
   By::AbstractVector{T}
   Bz::AbstractVector{T}
   Bmag::AbstractVector{T}
end

"""
	plotBSteady(flyby=8; file, dirobs, npict)

Plot quasi-steady state magnetic field comparison between simulation results from `file`
end observation from `dirobs`. Second norm is calculated to estimate the difference.
# Arguments
- `flyby::Int`: flyby index in [1,2,7,8,28,29].
"""
function plotBSteady(flyby=8; file="box_var_4_n00080000.out", dirobs=".", npict=1)

   # Read observation data
   fileobs = joinpath(dirobs, "Galileo_G"*string(flyby)*"_flyby_MAG.dat")
   t, df = getBObs(fileobs)
   BobsStrength = hypot.(df.Bx, df.By, df.Bz)

   # Read simulation data
   dirsim = dirname(file)
   if isempty(dirsim)
      data = readdata(file; npict)
   else
      data = readdata(basename(file); dir=dirsim, npict)
   end

   # Interpolate simulation data to observation data
   nx, nw = data.head.nx, data.head.nw

   x = data.x[:,1,1,1]
   y = data.x[1,:,1,2]
   z = data.x[1,1,:,3]

   w = data.w
   # Find the correct index of variables
   bx_ = findfirst(x->x=="bx", data.head.wnames)
   by_ = findfirst(x->x=="by", data.head.wnames)
   bz_ = findfirst(x->x=="bz", data.head.wnames)
   Bx = @view data.w[:,:,:,bx_]
   By = @view data.w[:,:,:,by_]
   Bz = @view data.w[:,:,:,bz_]

   # Gridded interpolation
   knots = (x, y, z)
   itpBx = interpolate(knots, Bx, Gridded(Linear()))
   BxSim = itpBx.(df.X,df.Y,df.Z)

   itpBy = interpolate(knots, By, Gridded(Linear()))
   BySim = itpBy.(df.X,df.Y,df.Z)

   itpBz = interpolate(knots, Bz, Gridded(Linear()))
   BzSim = itpBz.(df.X,df.Y,df.Z)

   BsimStrength = hypot.(BxSim, BySim, BzSim)

   # Visualization
   figure(figsize=(12,5))
   formatter = matplotlib.dates.DateFormatter("%H:%M")
   ax1 = subplot(411)
   plot(t, df.Bx, label=L"B_x")
   plot(t, BxSim)
   ax1.xaxis.set_major_formatter(formatter)
   legend()
   ax2 = subplot(412)
   plot(t, df.By, label=L"B_y")
   plot(t, BySim)
   ax2.xaxis.set_major_formatter(formatter)
   legend()
   ax3 = subplot(413)
   plot(t, df.Bz, label=L"B_z")
   plot(t, BzSim)
   ax3.xaxis.set_major_formatter(formatter)
   legend()
   ax4 = subplot(414)
   plot(t, BobsStrength, label=L"B")
   plot(t, BsimStrength)
   ax4.xaxis.set_major_formatter(formatter)
   legend()
   tight_layout()
   nothing
end

"""
    interp_simulation(file, dir, t, df, tStart, firstpict)

Return interpolated time-accurate B field stored in `file` from snapshot `firstpict` onto
timestamps `t` at satellite locations in `df`.
`tStart` is the starting DateTime for the synthetic satellite.
   
This is designed to handle any number of simulation snapshots. If the simulation time
interval is smaller than the observation's (which is usually the case), than the magnetic
field along the trajectory before time `tstart` is obtained from the first single snapshot
`firstpict`, and the field after the end time is obtained from the last single snapshot
ipict=npict. In this way we avoid the artificial jump due to periodic interpolation.
"""
function interp_simulation(file, dir, t, df, tStart, firstpict)

   data = readdata(file; dir)

   # Interpolate simulation data to observation data
   nx, nw = data.head.nx, data.head.nw
   npict = data.list.npictinfiles

   firstpict > npict && @error "firstpict out of range!"

   npict += 1 - firstpict

   kStart = findfirst(tStart .- t .< Second(0))

   tEnd  = tStart + Second(npict-1)
   kEnd = findfirst(tEnd .- t .< Second(0))

   nSim = kStart + (length(t) - kEnd + 1) + npict - 2
   BxSim = Vector{Float32}(undef, nSim)
   BySim = Vector{Float32}(undef, nSim)
   BzSim = Vector{Float32}(undef, nSim)
   tSim = Vector{DateTime}(undef, nSim)
   tSim[1:kStart] = t[1:kStart]

   # This is an approximation
   tSim[kStart:kStart+npict-1] .= timelinspace(tStart, tEnd, n=npict)
   tSim[kStart+npict-1:end] = t[kEnd:end]

   # Extract B field along the traj. before the starting time in one snapshot
   data = readdata(file; dir, npict=firstpict)

   # Interpolate simulation data to observation data
   x = data.x[:,1,1,1]
   y = data.x[1,:,1,2]
   z = data.x[1,1,:,3]

   bx_ = findfirst(x->x=="bx", data.head.wnames)
   by_ = findfirst(x->x=="by", data.head.wnames)
   bz_ = findfirst(x->x=="bz", data.head.wnames)
   Bx = @view data.w[:,:,:,bx_]
   By = @view data.w[:,:,:,by_]
   Bz = @view data.w[:,:,:,bz_]

   # Gridded interpolation
   knots = (x, y, z)
   itpBx = interpolate(knots, Bx, Gridded(Linear()))
   BxSim[1:kStart] = itpBx.(df.X[1:kStart],df.Y[1:kStart],df.Z[1:kStart])

   itpBy = interpolate(knots, By, Gridded(Linear()))
   BySim[1:kStart] = itpBy.(df.X[1:kStart],df.Y[1:kStart],df.Z[1:kStart])

   itpBz = interpolate(knots, Bz, Gridded(Linear()))
   BzSim[1:kStart] = itpBz.(df.X[1:kStart],df.Y[1:kStart],df.Z[1:kStart])

   # Extract the magnetic field along the trajectory from multiple snapshot
   for ipict = 1:npict-1
      println("ipict=$(round(ipict,digits=2))")
      data = readdata(file; dir, npict=firstpict+ipict)
      # Interpolate simulation data to observation data
      x = data.x[:,1,1,1]
      y = data.x[1,:,1,2]
      z = data.x[1,1,:,3]

      Bx = @view data.w[:,:,:,bx_]
      By = @view data.w[:,:,:,by_]
      Bz = @view data.w[:,:,:,bz_]

      k = argmin( abs.(tSim[kStart+ipict] .- t) )

      itpBx = interpolate(knots, Bx, Gridded(Linear()))
      BxSim[kStart+ipict] = itpBx(df.X[k], df.Y[k], df.Z[k])

      itpBy = interpolate(knots, By, Gridded(Linear()))
      BySim[kStart+ipict] = itpBy(df.X[k], df.Y[k], df.Z[k])

      itpBz = interpolate(knots, Bz, Gridded(Linear()))
      BzSim[kStart+ipict] = itpBz(df.X[k], df.Y[k], df.Z[k])
   end

   BxSim[kStart+npict:end] = itpBx.(df.X[kEnd+1:end], df.Y[kEnd+1:end], df.Z[kEnd+1:end])
   BySim[kStart+npict:end] = itpBy.(df.X[kEnd+1:end], df.Y[kEnd+1:end], df.Z[kEnd+1:end])
   BzSim[kStart+npict:end] = itpBz.(df.X[kEnd+1:end], df.Y[kEnd+1:end], df.Z[kEnd+1:end])

   BsimStrength = hypot.(BxSim, BySim, BzSim)

   Bsim = BField(tSim, BxSim, BySim, BzSim, BsimStrength)
end


function timelinspace(d1::DateTime, d2::DateTime; n=2)
   if n == 1
      return d1
   end
   Δ = d2 - d1
   T = typeof(Δ)
   δ = T(round(Int, Dates.value(Δ)/(n - 1)))
   d2 = d1 + δ*(n - 1)
   return d1:δ:d2
end

"Read observation data from `file` and return the timestamps and data."
function getBObs(file)

   df = CSV.File(file; header=2, delim=" ", ignorerepeated=true)
   t = DateTime.(df.yr, df.month, df.day, df.hr, df.min, floor.(Int, df.sec),
      floor.(Int, 1e3 .* (df.sec - floor.(df.sec))) )

   t, df
end

##
#plotBTimeAccurate("box_Hall_B_1200.outs";
#   ObsDir="/Users/hyzhou/Ganymede/GalileoData/galileomagdata",
#   SimDir="/Users/hyzhou/Ganymede/Hall_AMR3/GM",)

#=
# Plot time-accurate magnetic field comparisons.
flyby = 8
filename = ["box_Hall_B_1200.outs", "box_PIC_B_1200.outs"]
dir = "/Users/hyzhou/Documents/Computer/Julia/BATSRUS/VisAnaJulia"
firstpict, lastpict = 1, 1

t, df = getBObs()
BobsStrength = hypot.(df.Bx, df.By, df.Bz)

# Select the starting time for synthetic satellite
if flyby == 8
   #tStart = DateTime(1997,5,7,15,48,7)
   tStart = DateTime(1997,5,7,15,45,40)# PIC
   tStartHall = DateTime(1997,5,7,15,45,48)# Hall
elseif flyby == 28
   tStart = DateTime(2000,5,20,10,0,12)
end

Hall = interp_simulation(filename[1], dir, t,df, tStartHall, firstpict)

PIC = interp_simulation(filename[2], dir, t,df, tStart, firstpict)

fig, ax = plt.subplots(4,1,figsize=(10.0,5.0))
plt.rc("font", family="serif", size=14)

# Plot observation data
ax[1].plot(t, df.Bx, "k", label="Obs")
ax[2].plot(t, df.By, "k", label="Obs")
ax[3].plot(t, df.Bz, "k", label="Obs")
ax[4].plot(t, BobsStrength, "k", label="Obs")

# Plot simulation data
ax[1].plot(Hall.t, Hall.Bx, label="Hall MHD")
ax[2].plot(Hall.t, Hall.By, label="Hall MHD")
ax[3].plot(Hall.t, Hall.Bz, label="Hall MHD")
ax[4].plot(Hall.t, Hall.Bmag, label="Hall MHD")

ax[1].plot(PIC.t, PIC.Bx, label="MHD-EPIC")
ax[2].plot(PIC.t, PIC.By, label="MHD-EPIC")
ax[3].plot(PIC.t, PIC.Bz, label="MHD-EPIC")
ax[4].plot(PIC.t, PIC.Bmag, label="MHD-EPIC")

ax[1].set_ylabel("Bx [nT]")
ax[2].set_ylabel("By [nT]")
ax[3].set_ylabel("Bz [nT]")
ax[4].set_ylabel("B [nT]")
suptitle("Galileo G8 Flyby Magnetic Field")
for (i,a) in enumerate(ax)
   a.set_xlim(t[1500],t[end-2700])
   if i != 4
      a.axes.xaxis.set_ticklabels([])
   else
      formatter = matplotlib.dates.DateFormatter("%H:%M")
      a.xaxis.set_major_formatter(formatter)
   end
   a.tick_params(which="both", direction="in", top=true, right=true)
   a.minorticks_on()
   i == 1 &&
   a.legend(loc="lower left", bbox_to_anchor=(0.0, 0.86), ncol=3, frameon=false)
end

plt.subplots_adjust(hspace=0)
fig.align_ylabels(ax)
=#
