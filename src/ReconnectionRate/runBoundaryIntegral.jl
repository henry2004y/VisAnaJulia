# Control script for calculating CPCP for multiple timesteps.
# Using Channels for parallel tasks.
#
# To run this on a server without x-server, you must switch to the headless
# version of Paraview.
#
# Hongyang Zhou, hyzhou@umich.edu 03/19/2020

## Import packages
using Distributed, Glob
@everywhere using DelimitedFiles, ParallelDataTransfer
if !isdefined(Main, :ReconnectionRate)
   @everywhere include("utility.jl")
   @everywhere using .ReconnectionRate
end

## Function and global definitions

# The first parameter of a channel is the Type (e.g.: Int) and the second one,
# the maximum number of elements allowed in the channel.
const jobs    = RemoteChannel(()->Channel{Int}(32))
const results = RemoteChannel(()->Channel{Tuple}(32))

@everywhere function process_file(i)
   file = files[i]
   @info "processing file $file ..."
   # Create Python script with the correct input filename for tracing
   txt = read(filePy_orig[1], String)
   tracePy = filePy_orig[1][1:end-3]*string(i)*".py"
   fileTrace = "streamline"*lpad(i,4,"0")*".txt"
   fileStream = joinpath(fileDir, fileTrace)
   fileBoundaryXYZ = "boundary"*lpad(i,4,"0")*".txt"

   if isfile(fileDir*fileBoundaryXYZ)
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

   # Step 3: find magnetosphere boundary
   x, y, z = findBoundary(fileStream, method=1, nSmooth=100, DoPlot=false)
   saveBoundary(x, y, z, joinpath(fileDir,fileBoundaryXYZ))
   rm(fileDir*fileTrace, force=true)

   # Create Python script with the correct input filename for interpolation
   txt = read(filePy_orig[2], String)
   interpPy = filePy_orig[2][1:end-3]*string(i)*".py"
   fileBoundary = "boundary_value"*lpad(i,4,"0")*".txt"
   open(interpPy, "w") do f
      txt = replace(txt, "test.vtu" => file)
      txt = replace(txt, "boundary.txt" => fileBoundaryXYZ)
      txt = replace(txt, "boundary_value.txt" => fileBoundary)
      write(f, txt)
   end

   # Step 4: interpolate onto the boundary points (in Paraview)
   run(`$run_pvbatch $interpPy`)
   rm(interpPy)
end

function write_CPCP(filenames, fileDir="./")

   ϕ = Vector{Float32}(undef,length(filenames))

   for (i,file) in enumerate(filenames)
      fileBoundary = "boundary_value"*lpad(i,4,"0")*".txt"
      # Step 5: 2D line integral
      file_integration = joinpath(fileDir, fileBoundary)
      ϕ[i] = integrate_along_boundary_hall(file_integration, false)

      open(joinpath(fileDir,"CPCP.txt"), "a") do io
         i_end  = findfirst("_n",file)[1] - 1
         second = parse(Int32, file[i_end-1:i_end])
         minute = parse(Int32, file[i_end-3:i_end-2])
         timestep = 60*minute + second
         writedlm(io, [timestep round(ϕ[i],digits=4)])
      end
   end

end

@everywhere function do_work(jobs, results) # Define work function everywhere.
   while true
      job_id = take!(jobs)
      exec_time = @elapsed process_file(job_id)
      put!(results, (job_id, exec_time, myid()))
   end
end

function make_jobs(n)
   for i in 1:n
      put!(jobs, i)
   end
end

## Main
if Sys.isapple() # Sys.MACHINE might be better
   run_pvbatch = "pvbatch"
else
   run_pvbatch = "run_pvbatch"
end
# Step 1: create seeding points (same for all timesteps)
@everywhere dataDir = "../../GM/"
@everywhere fileDir = "data/"

!isdir(fileDir) && mkdir(fileDir)

createSeeds("seeds.txt", fileDir)
# Loop over multiple snapshots
files = basename.(glob("cut*.vtu", dataDir))

sendto(workers(), files=files)

@everywhere filePy_orig = ["streamtrace_save.py", "Interpolate_On_Boundary_batch.py"]

nFile = length(files)

@async make_jobs(nFile) # Feed the jobs channel with "n" jobs.

for p in workers() # Start tasks on the workers to process requests in parallel.
   @async remote_do(do_work, p, jobs, results) # Similar to remotecall.
end

@elapsed while nFile > 0 # Print out results.
   job_id, exec_time, wid = take!(results)
   println("$job_id finished $(round(exec_time; digits=2)) seconds "*"
   on worker $wid")
   global nFile -= 1
end

write_CPCP(files, fileDir)
