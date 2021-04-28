"""
Reconnection rate postprocessing related to ParaView.

Based on the tests, using concave (1) or convex (2) hull will effect the
potential drop calculation. Method 1 will at most give a value 10% smaller than
method 2.

In the near future, I want to call ParaView directly in Julia.
Hongyang Zhou, hyzhou@umich.edu 11/11/2019
modified on 03/16/2020
"""
module ReconnectionRate

using DelimitedFiles, NumericalIntegration, LinearAlgebra, Random, PyPlot
using ConcaveHull, LazySets
#using MATLAB

export Rg, e
export findBoundary, saveBoundary, createSeeds, integrate_along_boundary,
   integrate_along_boundary_hall


const Rg = 2634000        # [m], radius of Ganymede
const e  = 1.60217662e-19 # [C], electron charge


"""
	findBoundary(filename; DoPlot=false, method=1, nSmooth=3)

Identify boundary points from Paraview polyline tracing output.
`method` 1 is the convex hull of the boundary points, 2 is the concave hull
computed with KNN, and 3 is the Matlab `boundary` function call.  `nSmooth` is
the number of point considered in KNN.
To use 2 or 3,
you need to import the corresponding packages.
"""
function findBoundary(filename; DoPlot=false, method=1, nSmooth=3)

   fieldline = readdlm(filename, ',', Float32, header=true)
   data = fieldline[1]

   # count the total number of field lines
   nlinetotal = 0
   start = Int32[]
   for i = 1:size(data,1)
      if data[i,4] == 0.0 # IntegrationTime starts at 0.0
         nlinetotal += 1
         push!(start, i)
      end
   end

   line = Array{Float32, 2}(undef, 6, nlinetotal)

   for i = 1:nlinetotal-1
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
   status = zeros(Int8, nlinetotal)

   # Cluster seed points based on their end points
   for i = 1:nlinetotal
      if line[4,i]^2 + line[5,i]^2 + line[6,i]^2 < 1.1^2
         status[i] = 1
      end
   end

   pts_closed = line[1:2, status .== 1]
   pts_open = line[1:2, status .== 0]

   if method == 1 # LazySets, convex hull
      x = pts_closed[1,:]
      y = pts_closed[2,:]
      v = [[x[i], y[i]] for i in 1:length(x)]
      hull = convex_hull(v) # LazySets
      xy = hcat(hull...)
      xB = [xy[1,:]; xy[1,1]]
      yB = [xy[2,:]; xy[2,1]]
   elseif method == 2 # ConcaveHull
      ϵ = Float32(1e-6)
      Random.seed!(1234)
      x = [i + ϵ * rand(Float32) for i in pts_closed[1,:]]
      y = [i + ϵ * rand(Float32) for i in pts_closed[2,:]]
      v = [[x[i], y[i]] for i in 1:length(x)]
      hull = concave_hull(v, nSmooth) # ConcaveHull
      xB = [p[1] for p in hull.vertices]
      yB = [p[2] for p in hull.vertices]
   elseif method ==3 # Matlab
      x = pts_closed[1,:]
      y = pts_closed[2,:]
      k, area = mxcall(:boundary, 2, x, y)
      k = Int.(k)
      xB, yB = x[k], y[k]
   end
   zB = fill(2.0f0, length(xB))

   if DoPlot
      fig, ax = subplots(1,2)
      ax[1].scatter(pts_closed[1,:], pts_closed[2,:])
      ax[1].scatter(pts_open[1,:], pts_open[2,:])
      ax[2].scatter(xB, yB)
   end

   return xB, yB, zB
end

function saveBoundary(x,y,z, filename="boundary.txt")
   open(filename, "w") do io
      write(io, "\"X\",\"Y\",\"Z\"\n")
      writedlm(io, [x y z], ',')
   end
end

"Create seeds for field line tracing."
function createSeeds(filename="seeds.txt", dir="./")

   xMin, xMax = -2.0, 4.3 # xrange
   yMin, yMax = -3.0, 3.1 # yrange
   z = 2.0 # z cut plane

   dx, dy = 1/32, 1/32 # resolutions

   xRange = xMin:dx:xMax
   yRange = yMin:dy:yMax

   X = [x for x in xRange, y in yRange if 3.5^2 > (x-1.2)^2 + (y-0.05)^2 > 2.7^2]
   Y = [y for x in xRange, y in yRange if 3.5^2 > (x-1.2)^2 + (y-0.05)^2 > 2.7^2]
   Z = fill(z, size(X))

   open(joinpath(dir,filename), "w") do io
      write(io, "\"X\",\"Y\",\"Z\"\n")
      writedlm(io, [X Y Z], ',')
   end
end


"Integrate U cross B along the 2D boundary line."
function integrate_along_boundary(filename, DoPlot=false)
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

"Integrate Ue cross B along the magnetopause boundary."
function integrate_along_boundary_hall(filename, DoPlot=false)
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

end
