# Utility functions for extracting timings from runlog files.
#
# Hongyang Zhou, hyzhou@umich.edu

using DataFrames, Printf, Query

function extractTestInfo(File::AbstractString, compiler::AbstractString="unknown")

   nMPI = 0
   nThread = 0
   df = DataFrame(Name=String[], Time=Float64[], Proc=Int32[], Thread=Int32[],
      Compiler=String[])

   f = open(File)
   lines = readlines(f)
   keylines = String[]
   IsComplete,iEnd,iSeq = false,0,1

   UseTimingTree_I = occursin.("TREE",lines)
   UseTimingTree = any(x->x==true,UseTimingTree_I)

   i,iParam = 1,Int8(-1) # initialization
   while true
      if occursin("version", lines[i])
         iParam = i
         break
      end
      i += 1
   end

   for i in reverse(eachindex(lines))
      iStart = 0
      if !IsComplete && occursin("report", lines[i])
         IsComplete = true; continue
      elseif IsComplete && startswith(lines[i],"-------")
         iEnd = i-1; continue
      elseif IsComplete && startswith(lines[i],"BATSRUS")
         iStart = i
         push!(keylines, lines[[iParam;iStart;collect(iStart+2:iEnd)]]...); break
      end
   end

   if !isempty(keylines)
      parameterline = keylines[1]
      timelines = keylines[2:end]
   else
      println("runlog file $File not complete! Check it again!")
      close(f); return
   end

   nProc = parse(Int,parameterline[42:49])
   nThread = parse(Int,parameterline[66:69])
   #println("nProc=",nProc," nThread=",nThread) # Debug

   UseTimingTree ? iSeq = 4 : iSeq = 2

   for timeline in timelines
      testlist = split(timeline)
      push!(df,[testlist[1],parse(Float64,testlist[iSeq]),nProc,nThread,compiler])
   end
   close(f)

   # Debug
   #@show df

   return df

end


function postprocess(DirIn::AbstractString=".")

   dfAll = DataFrame(Name=String[], Time=Float64[], Proc=Int32[], Thread=Int32[],
      Compiler=String[])

   # Extract run time parameters and timings from runlog files
   for (root, dirs, files) in walkdir(DirIn)
      println("Directories in $root")
      #for dir in dirs
      #   println(joinpath(root, dir)) # path to directories
      #end
      println("Files in $root")
      for file in filter(x -> startswith(x,"runlog"),files)
         println(joinpath(root, file)) # path to files
         m = match(r"(intel|gnu|cray|pgi)",root)
         df = extractTestInfo(joinpath(root,file),m.match)
         try
            append!(dfAll,df)
         catch
            println("No test info from file ",file)
         end
      end
   end


   table = @from p in dfAll begin
      @where p.Name == "BATSRUS" #&& p.Proc*p.Thread == 256
      @select {p.Name, p.Time, p.Proc, p.Thread, p.Compiler}
      @collect DataFrame
   end

   return table
end
