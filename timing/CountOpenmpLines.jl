# Count the OpenMP changes to the BATS-R-US source code.
#
# Hongyang Zhou, hyzhou@umich.edu 05/09/2019

using DataFrames, Printf, Query

function extractOpenmpInfo(File::AbstractString; Dir::AbstractString=".",
   Debug::Bool=false)

nLine,nOpenmpLine = 0,0 # Initialize counters
for line in eachline(joinpath(Dir,File))
   if occursin(r"^\s*!\$omp",line) nOpenmpLine += 1 end
   # Exclude comment/empty lines
   if !(occursin(r"^\s*![^$]",line) || occursin(r"^\s*$",line))
      nLine += 1
   end
end

if Debug
   println("Searching in $File...")
   println("Total lines:",nLine)
   println("Total OpenMP lines:",nOpenmpLine)
end

df = DataFrame(Name=File, TotalLine=nLine, OpenmpLine=nOpenmpLine)

return df

end


function postprocess(DirIn::AbstractString=".")

dfFiles = DataFrame(Name=String[], TotalLine=Int32[], OpenmpLine=Int32[])

# Extract run time parameters and timings from runlog files
for (root, dirs, files) in walkdir(DirIn)
   !occursin(r"\/src\w*$",root) && continue

   println("Files in $root")
   for file in filter(x -> endswith(x,"f90"),files)
      df = extractOpenmpInfo(file,Dir=root)
      try
         append!(dfFiles,df)
      catch
         println("No info from file ",file)
      end
   end
end

return dfFiles

end

dfFiles = postprocess("/Users/hyzhou/SWMF/test/BATSRUS/srcBATL")
display(describe(dfFiles))
#showall(dfFiles)
ratio = sum(dfFiles[:OpenmpLine]) / sum(dfFiles[:TotalLine])
@printf("OpenMP ratio = %4.2f%%",ratio*100)
