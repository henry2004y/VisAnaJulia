# Scripts for field line postprocessing.
#
# Create seeds for field line tracing
# Identify boundary points from Paraview polyline tracing output
# Integrate electric fields along the magnetopause boundary
#
# In the near future, I want to call Paraview directly in Julia.
#
# Hongyang Zhou, hyzhou@umich.edu 11/11/2019

using DelimitedFiles, PyPlot, NumericalIntegration, LinearAlgebra
using LazySets
#using MATLAB
#using ConcaveHull

const Rg = 2634000        # [m], radius of Ganymede
const e  = 1.60217662e-19 # [C], electron charge

function findBoundary(filename, DoPlot::Bool=false)

   fieldline = readdlm(filename, ',', header=true)
   data = fieldline[1]

   # count the total number of field lines
   ilinetotal = 0
   start = []
   for i=1:size(data,1)
      if data[i,4] == 0.0
         ilinetotal += 1
         push!(start, i)
      end
   end

   line = Array{Float64, 2}(undef, 6, ilinetotal)
   iline = 1

   for i=1:ilinetotal-1
      iStart, iEnd = start[i], start[i+1]-1
      line[1,i] = data[iStart,5]
      line[2,i] = data[iStart,6]
      line[3,i] = data[iStart,7]
      line[4,i] = data[iEnd,5]
      line[5,i] = data[iEnd,6]
      line[6,i] = data[iEnd,7]
   end

   line[1,end] = data[start[end],5]
   line[2,end] = data[start[end],6]
   line[3,end] = data[start[end],7]
   line[4,end] = data[end,5]
   line[5,end] = data[end,6]
   line[6,end] = data[end,7]

   # status
   # (0: open, 1: closed along B, 2: closed along -B 3: fully closed)
   # (-1: cells inside body, -2: loop ray within block, -3: strange)
   status = fill(0, ilinetotal)

   # Cluster seed points based on their end points
   for i=1:ilinetotal
      if line[4,i]^2 + line[5,i]^2 + line[6,i]^2 < 1.1^2
         status[i] = 1
         #elseif abs(line[4,i]) ≥ 2.9 || abs(line[5,i]) ≥ 2.9 || abs(line[6,i]) ≥ 2.9
         #   status[i] = 0
      end
   end

   pts_closed = line[1:3, status .== 1]
   pts_open = line[1:3, status .== 0]

   if DoPlot
      scatter3D(pts_closed[1,:], pts_closed[2,:], pts_closed[3,:])
      scatter3D(pts_open[1,:], pts_open[2,:], pts_open[3,:])
   end

   x = pts_closed[1,:]
   y = pts_closed[2,:]
   z = pts_closed[3,:]

   # Matlab
   #=
   k, area = mxcall(:boundary, 2, x, y)
   k = Int.(k)
   figure()
   plot(x[k], y[k])
   return x[k], y[k], z[k]
   =#

   # LazySets
   v = [[x[i], y[i]] for i in 1:length(x)]
   hull = convex_hull(v) # LazySets
   xy = hcat(hull...)
   x = [xy[1,:]; xy[1,1]]
   y = [xy[2,:]; xy[2,1]]
   z = fill(2.0, length(x))
   return x, y, z

   # ConcaveHull
   #=
   v = [[x[i], y[i]] for i in 1:length(x)]
   hull = concave_hull(v,200) # ConcaveHull
   x = [p[1] for p in hull.vertices]
   y = [p[2] for p in hull.vertices]
   z = fill(2.0, length(x))
   return x,y,z
   =#
end

function saveBoundary(x,y,z, filename::String="boundary.txt")
   open(filename, "w") do io
      write(io, "\"X\",\"Y\",\"Z\"\n")
      writedlm(io, [x y z], ',')
   end
end

function createSeeds(filename::String="seeds.txt", dir::String="./", Verbose::Bool=false)

   xMin, xMax = -2.0, 4.3
   yMin, yMax = -3.0, 3.0
   z = 2.0

   dx, dy = 5e-2, 5e-2

   xRange = xMin:dx:xMax
   yRange = yMin:dy:yMax

   X = [x for x in xRange, y in yRange if (x-1)^2 + y^2 > 2.4^2 && (x-1)^2 + y^2 < 3.5^2]
   Y = [y for x in xRange, y in yRange if (x-1)^2 + y^2 > 2.4^2 && (x-1)^2 + y^2 < 3.5^2]
   Z = fill(z, size(X))

   open(dir*filename, "w") do io
      write(io, "\"X\",\"Y\",\"Z\"\n")
      writedlm(io, [X Y Z], ',')
   end

   Verbose && println("Seeds created.")
end


"""
   integrate_along_boundary(filename)
Integrate variables along the 2D boundary line.

May add Hall velocity later!
"""
function integrate_along_boundary(filename, DoPlot::Bool=false)
   bc = readdlm(filename, ',', header=true)
   data = bc[1]

   x = @view data[:,7]
   y = @view data[:,8]

   # Get E field: E = -uxB
   Ex = Vector{Float64}(undef,size(data,1))
   Ey = Vector{Float64}(undef,size(data,1))
   Ez = Vector{Float64}(undef,size(data,1))

   for i = 1:size(data,1)
      Ex[i], Ey[i], Ez[i] = -cross(data[i,1:3], data[i,4:6])
   end

   ϕx = cumul_integrate(x, Ex)
   ϕy = cumul_integrate(y, Ey)
   Coef = 1e3 * 1e-9 * Rg * 1e-3 # [km/s]*[nT]*[Rg] --> [kV]
   ϕ  = (ϕx .+ ϕy) .* Coef # [kV]

   if DoPlot
      figure()
      plot(ϕ)
   end

   maximum(ϕ) - minimum(ϕ)
end


function integrate_along_boundary_hall(filename, DoPlot::Bool=false)
   bc = readdlm(filename, ',', header=true)
   data = bc[1]

   ρ = @view data[:,1]
   ux= @view data[:,2]
   uy= @view data[:,3]
   uz= @view data[:,4]
   bx= @view data[:,5]
   by= @view data[:,6]
   bz= @view data[:,7]
   jx= @view data[:,8]
   jy= @view data[:,9]
   jz= @view data[:,10]
   x = @view data[:,11]
   y = @view data[:,12]

   # Get E field: E = -uxB
   Ex = Vector{Float64}(undef,size(data,1))
   Ey = Vector{Float64}(undef,size(data,1))
   Ez = Vector{Float64}(undef,size(data,1))

   # Calculate electron bulk velocity in Hall MHD
   # Transform the hall velocity into [km/s]
   Coef = 1e-6 /e * 1e-3 * 1e-6
   uex = @. ux - jx/ρ * Coef
   uey = @. uy - jy/ρ * Coef
   uez = @. uz - jz/ρ * Coef

   for i = 1:size(data,1)
      Ex[i], Ey[i], Ez[i] =
         -cross([uex[i], uey[i], uez[i]], [bx[i], by[i], bz[i]])
   end

   ϕx = cumul_integrate(x, Ex)
   ϕy = cumul_integrate(y, Ey)
   Coef = 1e3 * 1e-9 * Rg * 1e-3 # [km/s]*[nT]*[Rg] --> [kV]
   ϕ  = (ϕx .+ ϕy) .* Coef # [kV]

   if DoPlot
      figure()
      plot(ϕ)
   end

   #maximum(ϕ) - minimum(ϕ) # This would return the CPCP on the night side
   maximum(ϕ) - minimum(ϕ) + ϕ[end] # Return the CPCP on the day side
end
