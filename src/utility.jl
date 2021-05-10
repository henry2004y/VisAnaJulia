# Utility functions for space data analysis.

using ImageFiltering, OffsetArrays, FFTW, Statistics, LinearAlgebra

export mag2db, sma, spectrum

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
    sma(x, n=5)

Simple moving box average of the vector data `x` with box length 'n'.
One-sided average on the left and right edge with replicate border.
"""
function sma(x::Vector, n::Int=5)
   iseven(n) && @warn "Even box length detected!"
   sma(x, ((n-1)÷2, (n-1)÷2))
end

"""
    spectrum(x, Fs)

Return the frequency and amplitude for the single-sided spectrum of vector `x`
given sample frequency of `Fs`.
"""
function spectrum(x::Vector, Fs)
   L = length(x)
   x̃ = fft(x)
   # Compute the two-sided spectrum amplitude P2.
   P₂ = @. abs(x̃/L)
   # Compute the single-sided spectrum amplitude P1.
   P₁ = P₂[1:floor(Int, L/2+1)]
   P₁[2:end-1] .= 2 .* P₁[2:end-1]

   f = Fs * (0:(L/2)) / L

   f, P₁
end