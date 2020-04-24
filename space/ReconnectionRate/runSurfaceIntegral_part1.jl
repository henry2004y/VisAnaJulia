# Cross Polar Cap Potential (CPCP) calculation at z=2.0 cut plane.
# Combination of surface integral and line integral.
# Part 1: parallel ParaView launching on clusters for stream tracing and point
# interpolation.
#
# Using Channels for parallel tasks.
#
# Hongyang Zhou, hyzhou@umich.edu 03/19/2020

using Distributed, Glob
@everywhere using DelimitedFiles, ParallelDataTransfer

# The first parameter of a channel is the Type (e.g.: Int) and the second one,
# the maximum number of elements allowed in the channel.
const jobs    = RemoteChannel(()->Channel{Int}(32))
const results = RemoteChannel(()->Channel{Tuple}(32))

"Create surface seeds for interpolation."
function createSurfaceSeeds(filename="seeds.txt", dir="./")
   xMin, xMax = -2.0, 1.3 # xrange
   yMin, yMax = -3.0, 3.1 # yrange
   z = 2.0 # z cut plane

   dx, dy = 1/32, 1/32 # resolutions

   xRange = xMin:dx:xMax
   yRange = yMin:dy:yMax

   X = [x for x in xRange, y in yRange if 3.5^2 > (x-1.2)^2 + (y-0.05)^2]
   Y = [y for x in xRange, y in yRange if 3.5^2 > (x-1.2)^2 + (y-0.05)^2]
   Z = fill(z, size(X))

   open(joinpath(dir,filename), "w") do io
      write(io, "\"X\",\"Y\",\"Z\"\n")
      writedlm(io, [X Y Z], ',')
   end
end

@everywhere function findStatus(filename, fileDir, outname)

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
      if line[4,i]^2 + line[5,i]^2 + line[6,i]^2 < 1.2^2 ||
         (line[4,i] > 0.0 && (-2.0 < line[5,i] < 2.0)) # no holes inside!
         status[i] = 1
      end
   end

   pts_closed = line[1:2, status .== 1]
   x = @view pts_closed[1,:]
   y = @view pts_closed[2,:]
   z = fill(2.0f0, length(x))

   open(fileDir*outname, "w") do io
      write(io, "\"X\",\"Y\",\"Z\"\n")
      writedlm(io, [x y z], ',')
   end
end

@everywhere function process_file(i)

   file = files[i]
   @info "processing file $file ..."
   # Create Python script with the correct input filename for tracing
   txt = read(filePy_orig[1], String)
   tracePy = filePy_orig[1][1:end-3]*string(i)*".py"
   fileTrace = "streamline"*lpad(i,4,"0")*".txt"
   fileStream = joinpath(fileDir, fileTrace)
   fileSurfaceXYZ = "surface"*lpad(i,4,"0")*".txt"
   interpPy = filePy_orig[2][1:end-3]*string(i)*".py"
   fileSurface = "surface_value"*lpad(i,4,"0")*".txt"

   if isfile(fileDir*fileSurface)
      return # skip the finished snapshots
   end

   open(tracePy, "w") do f
      txt = replace(txt, "test.vtu"=>file)
      txt = replace(txt, "streamline.txt"=>fileTrace)
      write(f, txt)
   end

   # Step 2: field line tracing (in Paraview)
   run(`$run_pvbatch $tracePy`)
   rm(tracePy, force=true)

   # Step 3: find open/closed status for each point
   findStatus(fileStream, fileDir, fileSurfaceXYZ)
   rm(fileDir*fileTrace, force=true) # hang if file not found with force=false

   # Create Python script with the correct input filename for interpolation
   txt = read(filePy_orig[2], String)
   open(interpPy, "w") do f
      txt = replace(txt, "test.vtu" => file)
      txt = replace(txt, "surface.txt" => fileSurfaceXYZ)
      txt = replace(txt, "surface_value.txt" => fileSurface)
      write(f, txt)
   end

   # Step 4: interpolate onto the boundary points (in Paraview)
   run(`$run_pvbatch $interpPy`)
   rm(interpPy, force=true)
end

@everywhere function do_work(jobs, results) # Define work function everywhere.
   while true
      job_id = take!(jobs)
      process_file(job_id)
      put!(results, (job_id, myid()))
   end
end

function make_jobs(n)
   for i in 1:n
      put!(jobs, i)
   end
end


## Main
if Sys.isapple() # Sys.MACHINE might be better
   @everywhere run_pvbatch = "pvbatch"
else
   @everywhere run_pvbatch = "run_pvbatch"
end
##
@everywhere dataDir = "../../GM/"
@everywhere fileDir = "data/"

if !isdir(fileDir)
   mkdir(fileDir)
end
# Step 1: create seeding points (same for all timesteps)
createSurfaceSeeds("seeds.txt", fileDir)
# Loop over multiple snapshots
files = basename.(glob("cut*.vtu", dataDir))

sendto(workers(), files=files)

@everywhere filePy_orig = ["streamtrace_save.py", "Interpolate_On_Surface_batch.py"]

nFile = length(files)

@async make_jobs(nFile) # Feed the jobs channel with "n" jobs.

# Start tasks on the workers to process requests in parallel.
for p in workers()
   @async remote_do(do_work, p, jobs, results) # Similar to remotecall.
end

@elapsed while nFile > 0 # Print out results.
   job_id, wid = take!(results)
   println("$job_id finished on worker $wid")
   global nFile -= 1
end

##
# Integration after the parallel work using ParaView is finished.
# See runSurfaceIntegral_part2.jl
