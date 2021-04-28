# Utility functions for space data analysis.

using ImageFiltering, OffsetArrays

export sma, read_satellite_data, mag2db

"Import satellite data. Require one line header presented."
function read_satellite_data(fname)
   f = readdlm(fname, ',', Float32, '\n'; header=true)

   header = f[2][1:end-1]
   data   = f[1]

   satelliteNo = unique(data[:,1]) # number of static satellites

   if minimum(satelliteNo) > 0.1
      @error "Some thing is wrong with the input data!"
   end

   header, data, satelliteNo
end

"Convert `x` from magnitude to decibels."
mag2db(x) = 10. * log10.(x)

"""
    sma(x, n=(5,5))

Simple moving average of `x` with `n[1]` leading and `n[2]` trailing points.
Default boundary is "replicate", meaning that the border pixels extend beyond
the boundaries.
"""
function sma(x, n::Tuple{Int,Int}=(5,5))
   len = sum(n) + 1
   kernel = OffsetArray(fill(1/len, len), -n[1]:n[2])
   imfilter(x, kernel)
end


"""
    sma(x, n=100)

Simple moving box average of the vector data `x` with box length 'n'.
One-sided average on the left and right edge with replicate border.
"""
function sma(x::Vector, n::Int=5)
   iseven(n) && @warn "Even box length detected!"
   sma(x, ((n-1)÷2, (n-1)÷2))
end


"""
    ema(x, n=5)

Return the exponentially moving box average of the vector data `x` with box
length 'n'.
Check https://github.com/JuliaQuant/MarketTechnicals.jl/blob/master/src/movingaverages.jl
"""
function ema(x::Vector, n::Int=5)
   @warn "To be implemented!"
end
