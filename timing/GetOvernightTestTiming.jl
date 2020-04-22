# Web Scraping in Julia
# Getting nightly test information from web.
#
#
# Hongyang Zhou, hyzhou@umich.edu 11/22/2018 initial version
# 05/08/2019 add timings for different tests; add parameter input file

# I need a way to check if the test correctly passed!
# The logic imposed is really fragile...

using HTTP, Gumbo, Dates, DelimitedFiles, ProgressMeter, TOML
using Plots
using Statistics
using Printf
pyplot() # Set the backend to Plotly
#plotly()

"""
   getFuncTimings(dates,[testsuite,testcase])

Obtain the timings from html files for SWMF tests.

# Arguments
- `dates::Vector{Date}`: date interval.
- `testsuite::String="mesh"`: name of machine/platform.
"""
function getFuncTimings(dates::Array{Date,1}; testsuite::String = "mesh",
   verbose::Bool=false)
   url_front = "http://herot.engin.umich.edu/~gtoth/SWMF_TEST_RESULTS/"
   url_end   = "/test_swmf.html"
   ndates = size(dates,1)
   urls = Array{String}(undef,ndates)
   timings = zeros(ndates)
   tTotal = Vector{Float64}(undef,ndates)

   # Obtain data from html
   @showprogress 1 "Fetching data..." for i=1:ndates
      urls[i] = url_front*Dates.format(dates[i],"yyyy/mm/dd/")*testsuite*url_end

      try
         r = HTTP.request("GET", urls[i]; retries=0,verbose=0)
         rawText = String(r.body)[end-2850:end]
         rawText = rawText[ collect(findfirst("\n000",rawText))[1]+1:end ]

         timeTable = split(rawText,"\n")
         for row in timeTable[1:36]
            cols = split(row)
            try
               timings[i] = parse(Float64,cols[3])
            catch
               verbose && println("Test",cols[1]," failed on ",dates[i])
               continue
            end
         end

         tTotal[i] = parse(Float64,split(timeTable[38])[3])
      catch
         continue
      end

   end

   # Clear all large/small values
   for i=1:length(tTotal)
      if tTotal[i] < 1e-5 || tTotal[i] > 1e4
         tTotal[i] = 0.0
      elseif isnan(tTotal[i])
         tTotal[i] = 0.0
      end
   end

   return timings,tTotal
end

"""
   getTimings(dates,[testsuite,testcase])

Obtain the timings from html log files for SWMF tests.

# Arguments
- `dates::Vector{Date}`: date interval.
- `testsuite::String="mesh"`: name of machine/platform.
"""
function getTimings(dates::Array{Date,1};
   testsuite::String="gfortran",testcase::String="shocktube",
   verbose::Bool=false)
# Find test_shocktube_run
# PostProc.pl: TIMINGS from runlog (init, run) 0.01 1.28
url_front = "http://herot.engin.umich.edu/~gtoth/SWMF_TEST_RESULTS/"
url_end   = "/test_swmf.log"
ndates = size(dates,1)
urls = Array{String}(undef,ndates)
timings = zeros(ndates)
nFail,nError = 0,0

# Obtain data from log
@showprogress 1 "Fetching data..." for i=1:ndates
   urls[i] = url_front*Dates.format(dates[i],"yyyy/mm/dd/")*testsuite*url_end

   try
      r = HTTP.request("GET", urls[i]; retries=0,verbose=0)
      rawText = String(r.body)
      start_  = collect(findfirst(testcase*"_run\n",rawText))[end]
      end_    = collect(findfirst(testcase*".diff\n",rawText))[1] + 100

      if end_ - start_ > 5000 # something went wrong if indexes diff to much
         timings[i] = 0.0
         nFail += 1
         continue
      end

      rawText = rawText[start_:end_]

      start_  = findfirst("(init, run)",rawText)[end]+7
      end_    = start_ + 3

      timings[i] = parse(Float64,rawText[start_:end_])

      end_ = collect(findfirst(testcase*".diff\n",rawText))[end]
      end_ = collect(findnext(testcase*".diff\n",rawText,end_))[1]

      start_ = end_ - 21
      end_ -= 19

      if DoCheckResult && parse(Int,rawText[start_:end_]) != 0
         timings[i] = 0.0
         nFail += 1
      end

   catch
      nError += 1
      continue
   end

end

# Remove outliers
tMean,tStd = mean(timings),std(timings)
for i=1:ndates
   if timings[i] > tMean + 5*tStd
      timings[i] = 0.0
   end
end

@printf("test pass rate = %4.2f%%, no data rate = %4.2f%%",
   (1.0-nFail/(ndates-nError))*100, nError/ndates*100)

return timings

end


paramIn = TOML.parsefile("PARAM.toml")
verbose   = paramIn["Param"]["verbose"]
testsuite = paramIn["Param"]["testsuite"]
testcase  = paramIn["Param"]["testcase"]
DoCheckResult = paramIn["Param"]["DoCheckResult"]

dateStart = Date(paramIn["Param"]["dateStart"])
dateEnd   = Date(paramIn["Param"]["dateEnd"])
dateGap   = Day(paramIn["Param"]["dateGap"])
dr = collect(dateStart:dateGap:dateEnd)


#timings, tTotal = getFuncTimings(dr)
#plot(dr,tAll[:,4,3], linewidth=2, label=labels, legend=:best)
#p1 = plot(dr,tTotal, linewidth=2,xtickfontsize=16, ytickfontsize=16,
#   title="Total Test Func Timing, mesh nagfor mac",
#   xlabel="date",ylabel="execution time [s]")

#=
timings, tTotal = getFuncTimings(dr,testsuite="gfortran",verbose=verbose)

p2 = plot(dr,tTotal, linewidth=2, marker=(:dot), markersize=1,
   xtickfontsize=16, ytickfontsize=16,
   title="Total Test Func Timing, gfortran mac",
   xlabel="date",ylabel="execution time [s]")
=#
if paramIn["Param"]["timeType"] == 1

   timings, tTotal = getFuncTimings(dr,testsuite="pleiades",verbose=verbose)

   p3 = plot(dr,tTotal, linewidth=2, marker=(:dot), markersize=1,
      label="FuncTest total", xtickfontsize=16, ytickfontsize=16,
      title="Total Test Func Timing, pleiades ifort",
      xlabel="date",ylabel="execution time [s]")
end


if paramIn["Param"]["timeType"] == 2

   timings = getTimings(dr,testsuite=testsuite,testcase=testcase)

   p4 = plot(dr,timings, linewidth=2, marker=(:dot), markersize=1,
      label=testcase, xtickfontsize=16, ytickfontsize=16,
      title="Nightly test timing, "*testsuite,
      xlabel="date",ylabel="execution time [s]")
end
