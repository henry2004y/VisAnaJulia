# Script for SWMF Galileo Flyby Simulation Comparison
#
# Hongyang Zhou, hyzhou@umich.edu 02/03/2020

using SWMF, CSV, Dates, PyPlot, Interpolations

#=
## Calculate differences & save to output file
Bdiff = Bobs - Bsim
# There`s a difference between 2 definitions of 2-norm.
norm2 = mean(abs(Bdiff(~isnan(Bdiff))).^2)^(1/2)
# Find the index and value of the best parameter set
[val, idx] = min(norm2)
=#


"""
	plotBSteady(flyby=8, DoSave=true)

Plot quasi-steady state magnetic field comparison. Second norm is calculated to
quantitatively estimate the difference.
`flyby` is an integer in [1,2,7,8,28,29].
`DoSave` decides
"""
function plotBSteady(flyby=8; filename="box_var_4_n00080000.out", DoSave=false,
   dir="/Users/hyzhou/Documents/research/Ganymede/Galileo")
   # Read observation data

   flybyfile = "Galileo_G"*string(flyby)*"_flyby_MAG.dat"

   df = CSV.File(joinpath(dir,flybyfile);header=2,delim=" ",ignorerepeated=true)
   t = DateTime.(df.yr, df.month, df.day, df.hr, df.min, floor.(Int, df.sec),
   	floor.(Int, 1e3 .* (df.sec - floor.(df.sec))) )
   BobsStrength = @. sqrt(df.Bx^2 + df.By^2 + df.Bz^2)

   # Read simulation data

   npict = 1 # Remember to change this for different runs!
   head, data = readdata(filename,
      dir="/Users/hyzhou/Documents/Computer/Julia/BATSRUS/VisAnaJulia",
      npict=npict)

   # Interpolate simulation data to observation data
   nx, nw = head[1][:nx], head[1][:nw]

   x = data[1].x[:,1,1,1]
   y = data[1].x[1,:,1,2]
   z = data[1].x[1,1,:,3]

   w = data[1].w
   # Find the correct index of variables
   bx_ = findfirst(x->x=="bx", head[1][:wnames])
   by_ = findfirst(x->x=="by", head[1][:wnames])
   bz_ = findfirst(x->x=="bz", head[1][:wnames])
   Bx = @view data[1].w[:,:,:,bx_]
   By = @view data[1].w[:,:,:,by_]
   Bz = @view data[1].w[:,:,:,bz_]

   # Gridded interpolation
   knots = (x, y, z)
   itpBx = interpolate(knots, Bx, Gridded(Linear()))
   BxSim = itpBx.(df.X,df.Y,df.Z)

   itpBy = interpolate(knots, By, Gridded(Linear()))
   BySim = itpBy.(df.X,df.Y,df.Z)

   itpBz = interpolate(knots, Bz, Gridded(Linear()))
   BzSim = itpBz.(df.X,df.Y,df.Z)

   BsimStrength = @. sqrt(BxSim^2 + BySim^2 + BzSim^2)

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

   if DoSave # Save the plot

   end

end

"""
	plotBTimeAccurate(filename, flyby=8, firstpict=1, lastpict=1;
   	DoSave=false, IsOneFile=false, ObsDir, SimDir)

Plot time-accurate magnetic field comparisons.
This function is able to handle any npict number of snapshots. If the simulation
time interval is smaller than observation (which is usually the case), than the
magnetic field along the trajectoryy before time_start is obtained from first
single snapshot ipict=1, and B after time_end is obtained from the last single
snapshot ipict=npict. In this way, I avoid the jump due to periodic
interpolation.
`flyby` and `time_start` need to be modified when switching between flybys.
"""
function plotBTimeAccurate(filename::String, flyby=8, firstpict=1, lastpict=1;
   DoSave=false, IsOneFile=false,
   ObsDir="/Users/hyzhou/Documents/research/Ganymede/Galileo",
   SimDir="/Users/hyzhou/Documents/Computer/Julia/BATSRUS/VisAnaJulia")

   # Read observation data
   Obsfile = "Galileo_G"*string(flyby)*"_flyby_MAG.dat"

   df = CSV.File(joinpath(ObsDir,Obsfile);header=2,delim=" ",ignorerepeated=true)
   t = DateTime.(df.yr, df.month, df.day, df.hr, df.min, floor.(Int, df.sec),
      floor.(Int, 1e3 .* (df.sec - floor.(df.sec))) )
   BobsStrength = @. sqrt(df.Bx^2 + df.By^2 + df.Bz^2)

   # Select the starting time for synthetic satellite
   if flyby == 8
      tStart = DateTime(1997,5,7,15,48,7)
   elseif flyby == 28
      tStart = DateTime(2000,5,20,10,0,12)
   end

   head, data, list = readdata(filename, dir=SimDir)

   # Interpolate simulation data to observation data
   nx, nw = head[1][:nx], head[1][:nw]
   npict = list[1].npictinfiles

   firstpict > npict && @error "firstpict out of range!"

   npict += 1 - firstpict

   kStart = findfirst(tStart .- t .< Second(0))

   tEnd  = tStart + Second(npict-1)
   kEnd = findfirst(tEnd .- t .< Second(0))

   nSim = kStart + (length(t) - kEnd + 1) + npict - 2
   BxSim = Vector{Float32}(undef,nSim)
   BySim = Vector{Float32}(undef,nSim)
   BzSim = Vector{Float32}(undef,nSim)
   tSim = Vector{DateTime}(undef,nSim)
   tSim[1:kStart] = t[1:kStart]

   # This is an approximation
   tSim[kStart:kStart+npict-1] .= timelinspace(tStart, tEnd, n=npict)
   tSim[kStart+npict-1:end] = t[kEnd:end]

   # Extract B field along the traj. before the starting time in one snapshot
   head, data = readdata(filename, dir=SimDir, npict=firstpict)

   # Interpolate simulation data to observation data
   x = data[1].x[:,1,1,1]
   y = data[1].x[1,:,1,2]
   z = data[1].x[1,1,:,3]

   bx_ = findfirst(x->x=="bx", head[1][:wnames])
   by_ = findfirst(x->x=="by", head[1][:wnames])
   bz_ = findfirst(x->x=="bz", head[1][:wnames])
   Bx = @view data[1].w[:,:,:,bx_]
   By = @view data[1].w[:,:,:,by_]
   Bz = @view data[1].w[:,:,:,bz_]

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
      head, data = readdata(filename, dir=SimDir, npict=firstpict+ipict)
      # Interpolate simulation data to observation data
      x = data[1].x[:,1,1,1]
      y = data[1].x[1,:,1,2]
      z = data[1].x[1,1,:,3]

      Bx = @view data[1].w[:,:,:,bx_]
      By = @view data[1].w[:,:,:,by_]
      Bz = @view data[1].w[:,:,:,bz_]

      k = argmin( abs.(tSim[kStart+ipict] .- t) )

      itpBx = interpolate(knots, Bx, Gridded(Linear()))
      BxSim[kStart+ipict] = itpBx(df.X[k],df.Y[k],df.Z[k])

      itpBy = interpolate(knots, By, Gridded(Linear()))
      BySim[kStart+ipict] = itpBy(df.X[k],df.Y[k],df.Z[k])

      itpBz = interpolate(knots, Bz, Gridded(Linear()))
      BzSim[kStart+ipict] = itpBz(df.X[k],df.Y[k],df.Z[k])
   end


   BxSim[kStart+npict:end] = itpBx.(df.X[kEnd+1:end],df.Y[kEnd+1:end],df.Z[kEnd+1:end])
   BySim[kStart+npict:end] = itpBy.(df.X[kEnd+1:end],df.Y[kEnd+1:end],df.Z[kEnd+1:end])
   BzSim[kStart+npict:end] = itpBz.(df.X[kEnd+1:end],df.Y[kEnd+1:end],df.Z[kEnd+1:end])

   BsimStrength = @. sqrt(BxSim^2 + BySim^2 + BzSim^2)

   # Visualization
   figure(figsize=(12,8))
   formatter = matplotlib.dates.DateFormatter("%H:%M")
   ax1 = subplot(411)
   plot(t, df.Bx, label=L"B_x")
   plot(tSim, BxSim)
   ax1.xaxis.set_major_formatter(formatter)
   legend()
   ax2 = subplot(412)
   plot(t, df.By, label=L"B_y")
   plot(tSim, BySim)
   ax2.xaxis.set_major_formatter(formatter)
   legend()
   ax3 = subplot(413)
   plot(t, df.Bz, label=L"B_z")
   plot(tSim, BzSim)
   ax3.xaxis.set_major_formatter(formatter)
   legend()
   ax4 = subplot(414)
   plot(t, BobsStrength, label=L"B")
   plot(tSim, BsimStrength)
   ax4.xaxis.set_major_formatter(formatter)
   legend()
   tight_layout()

   # Save plots
   if DoSave

   end
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

plotBTimeAccurate(;ObsDir="/Users/hyzhou/Ganymede/GalileoData/galileomagdata",
   SimDir="/Users/hyzhou/Ganymede/Hall_AMR3/GM", filename="box_Hall_B_1200.outs")
